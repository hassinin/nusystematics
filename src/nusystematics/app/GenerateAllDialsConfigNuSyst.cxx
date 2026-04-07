// GenerateAllDialsConfigNuSyst.cxx
// Builds a "kitchen sink" generated_systematic_provider_configuration fhicl
// containing every available dial:
//   - all GENIE Reweight dials enumerated from genie::rew::GSyst_t at runtime
//     (no hardcoded dial names — discovered via GSyst::AsString)
//   - all standalone nusystematics provider tool configs found in --fcl-dir,
//     each loaded and validated; ones that fail to instantiate are reported
//     and skipped.
//
// Modes (selected via --mode): genierw, providers, all  (default: all)

#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/utility/make_instance.hh"
#include "nusystematics/utility/response_helper.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/parse.h"

#include "RwFramework/GSyst.h"

#include <algorithm>
#include <dirent.h>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace systtools;
using namespace nusyst;

// ===== CLI =================================================================
namespace cliopts {
std::string mode = "all"; // genierw | providers | all
std::string outputfile;
std::string fcl_dir;
std::string variation_descriptor = "[-3,-2,-1,0,1,2,3]";
bool include_zexp = false;
std::vector<std::string> skip_providers;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Generates a 'generated_systematic_provider_configuration' fhicl\n"
    "  containing every available dial.\n\n"
    "  Optional:\n"
    "    --mode <m>             genierw | providers | all (default: all)\n"
    "    -o <output.fcl>        Output file (default: stdout)\n"
    "    --fcl-dir <dir>        Tool-config fcl directory\n"
    "                           (default: $NUSYST/fcl)\n"
    "    --variation-descriptor \"[-3,-2,-1,0,1,2,3]\"\n"
    "                           Variation descriptor used for every GENIE RW\n"
    "                           dial (default shown).\n"
    "    --include-zexp         Use Z-expansion CCQE form factor dials\n"
    "                           instead of the dipole CCQE dials\n"
    "                           (mutually exclusive in nusystematics).\n"
    "    --skip <a,b,...>       Comma-separated tool-config filenames to\n"
    "                           skip in providers mode\n"
    "                           (e.g. CCQERPA.ToolConfig.fcl).\n"
    "    -?|--help              Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "--mode") cliopts::mode = argv[++opt];
    else if (s == "-o") cliopts::outputfile = argv[++opt];
    else if (s == "--fcl-dir") cliopts::fcl_dir = argv[++opt];
    else if (s == "--variation-descriptor") cliopts::variation_descriptor = argv[++opt];
    else if (s == "--include-zexp") cliopts::include_zexp = true;
    else if (s == "--skip") {
      std::string tok; std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::skip_providers.push_back(tok);
    } else {
      std::cout << "[ERROR]: Unknown option: " << s << std::endl;
      SayUsage(argv); exit(1);
    }
    opt++;
  }
  if (cliopts::mode != "genierw" && cliopts::mode != "providers" && cliopts::mode != "all") {
    std::cout << "[ERROR]: --mode must be one of: genierw, providers, all" << std::endl;
    exit(1);
  }
}

// ===== Helpers =============================================================
bool IsValidGENIEDialName(const std::string &name) {
  if (name.empty() || name == "-") return false;
  if (name.find("INVALID") != std::string::npos) return false;
  return true;
}

// Mutual-exclusion / redundancy filter for the GENIEReWeight provider.
// Returns true if the dial should be included in the default tool config.
//
// nusystematics' GENIEReWeight has constraints encoded in
// GENIEReWeightParamConfig.cc that we have to respect:
//
//   1. Dipole vs Z-expansion CCQE axial form factor dials are mutually
//      exclusive. Default: dipole. With --include-zexp: Z-expansion.
//
//   2. shape+norm Ma/Mv dials and shape-only twins reweight the same physics.
//      We pick the shape-only twins so that the corresponding Norm* dials
//      (NormCCQE, NormCCRES, NormNCRES) can also be included independently.
//      The {channel}IsShapeOnly tool_option flags are set to true accordingly.
//
//   This gives the maximum number of independent reweighting knobs that
//   nusystematics' GENIEReWeight will accept in a single instance.
bool DialIsActive(const std::string &name) {
  static const std::set<std::string> zexp_axial = {
    "ZNormCCQE", "ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE",
    "AxFFCCQEshape"
  };
  // The shape+norm dials we drop in favour of the shape-only twins + Norm*
  static const std::set<std::string> shape_norm_combos_to_drop = {
    "MaCCQE", "MaCCRES", "MvCCRES", "MaNCRES", "MvNCRES",
    "AhtBY", "BhtBY", "CV1uBY", "CV2uBY", "E0CCQE"
  };

  if (cliopts::include_zexp) {
    // In Z-exp mode, drop the dipole CCQE axial dials (and the shape+norm
    // combos which contain MaCCQE).
    if (shape_norm_combos_to_drop.count(name)) return false;
    // Also drop NormCCQE-related dials in Z-exp mode (Z-exp has its own norm)
    if (name == "NormCCQE" || name == "NormCCQEenu") return false;
    if (name == "MaCCQEshape" || name == "VecFFCCQEshape") return false;
  } else {
    // Dipole mode: drop Z-expansion dials
    if (zexp_axial.count(name)) return false;
    // Drop shape+norm Ma/Mv twins; we use the *shape variants + Norm*
    if (shape_norm_combos_to_drop.count(name)) return false;
  }
  return true;
}

// Build the GENIEReWeight tool-config ParameterSet by enumerating GSyst_t.
// Returns the parameter set (not yet wrapped in a top-level config).
fhicl::ParameterSet BuildGENIEReWeightToolConfig(std::vector<std::string> &included_dials,
                                                  std::vector<std::string> &skipped_dials) {
  fhicl::ParameterSet ps;
  ps.put<std::string>("tool_type", "GENIEReWeight");
  ps.put<std::string>("instance_name", "All");
  // Required to allow splineable parameters individually
  ps.put<bool>("ignore_parameter_dependence", true);
  ps.put<std::string>("genie_tune_name", "${GENIE_XSEC_TUNE}");
  // Use shape-only mode for the channels with separate Norm* dials so we can
  // include both the shape twin (Ma*shape) and the Norm* dial independently.
  if (!cliopts::include_zexp) {
    ps.put<bool>("MaCCQEIsShapeOnly", true);
  }
  ps.put<bool>("CCRESIsShapeOnly", true);
  ps.put<bool>("NCRESIsShapeOnly", true);
  ps.put<bool>("DISBYIsShapeOnly", true);

  for (int i = static_cast<int>(genie::rew::kNullSystematic) + 1;
       i < static_cast<int>(genie::rew::kNTwkDials); ++i) {
    genie::rew::GSyst_t s = static_cast<genie::rew::GSyst_t>(i);
    std::string name = genie::rew::GSyst::AsString(s);
    if (!IsValidGENIEDialName(name)) {
      skipped_dials.push_back(name + " (invalid name)");
      continue;
    }
    if (!DialIsActive(name)) {
      skipped_dials.push_back(name + " (excluded by FF mode)");
      continue;
    }
    ps.put<double>(name + "_central_value", 0.0);
    ps.put<std::string>(name + "_variation_descriptor", cliopts::variation_descriptor);
    included_dials.push_back(name);
  }

  return ps;
}

// Detect a tool-config-style fhicl: must contain a top-level "syst_providers" key.
bool IsToolConfigFile(const std::string &path) {
  std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
  try {
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(path, *fm);
    return ps.has_key("syst_providers");
  } catch (...) {
    return false;
  }
}

// List *.fcl files in a directory (non-recursive).
std::vector<std::string> ListFcls(const std::string &dir) {
  std::vector<std::string> out;
  DIR *d = opendir(dir.c_str());
  if (!d) return out;
  struct dirent *ent;
  while ((ent = readdir(d)) != nullptr) {
    std::string name = ent->d_name;
    if (name.size() < 5) continue;
    if (name.substr(name.size() - 4) != ".fcl") continue;
    out.push_back(dir + "/" + name);
  }
  closedir(d);
  std::sort(out.begin(), out.end());
  return out;
}

// Verify that a provider can be re-instantiated from its own parameter
// headers via response_helper-style loading. This catches providers (like
// CCQERPA) that successfully build their metadata via ConfigureFromToolConfig
// but then fail in SetupResponseCalculator when actually used downstream.
bool VerifyProviderLoad(IGENIESystProvider_tool *prov, std::string &err_out) {
  try {
    fhicl::ParameterSet inner;
    inner.put(prov->GetFullyQualifiedName(), prov->GetParameterHeadersDocument());
    inner.put<std::vector<std::string>>("syst_providers",
                                         {prov->GetFullyQualifiedName()});
    nusyst::response_helper rh;
    rh.LoadProvidersAndHeaders(inner);
    return true;
  } catch (std::exception &e) {
    err_out = e.what();
    return false;
  } catch (...) {
    err_out = "unknown exception";
    return false;
  }
}

// Try to instantiate providers from a tool-config fhicl. Returns the
// configured providers (empty on failure) and writes any error to err_out.
// Also verifies each provider can be loaded from its own parameter headers.
// Uses syst_param_id_offset to make sure paramIds remain unique across calls.
std::vector<std::unique_ptr<IGENIESystProvider_tool>>
TryLoadToolConfig(const std::string &path, paramId_t &syst_param_id_offset,
                  std::string &err_out) {
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
  try {
    std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(path, *fm);
    tools = systtools::ConfigureISystProvidersFromToolConfig<IGENIESystProvider_tool>(
        ps, nusyst::make_instance, "syst_providers", syst_param_id_offset);
  } catch (std::exception &e) {
    err_out = e.what();
    return tools;
  } catch (...) {
    err_out = "unknown exception";
    return tools;
  }

  // Verification pass: try to load each provider from its own parameter
  // headers. Drop providers that fail this stricter test.
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> verified;
  for (auto &tool : tools) {
    std::string verr;
    if (VerifyProviderLoad(tool.get(), verr)) {
      syst_param_id_offset += tool->GetSystMetaData().size();
      verified.push_back(std::move(tool));
    } else {
      err_out = "verify failed: " + verr;
      return std::vector<std::unique_ptr<IGENIESystProvider_tool>>{};
    }
  }
  return verified;
}

// Try to instantiate the GENIEReWeight provider from an in-memory tool-config.
// Also verifies it can round-trip through response_helper-style loading.
std::vector<std::unique_ptr<IGENIESystProvider_tool>>
TryLoadInMemoryGENIERW(const fhicl::ParameterSet &tool_ps,
                        paramId_t &syst_param_id_offset, std::string &err_out) {
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
  try {
    fhicl::ParameterSet wrapper;
    wrapper.put("GENIEReWeight_All", tool_ps);
    wrapper.put<std::vector<std::string>>("syst_providers", {"GENIEReWeight_All"});
    tools = systtools::ConfigureISystProvidersFromToolConfig<IGENIESystProvider_tool>(
        wrapper, nusyst::make_instance, "syst_providers", syst_param_id_offset);
  } catch (std::exception &e) {
    err_out = e.what();
    return tools;
  } catch (...) {
    err_out = "unknown exception";
    return tools;
  }

  std::vector<std::unique_ptr<IGENIESystProvider_tool>> verified;
  for (auto &tool : tools) {
    std::string verr;
    if (VerifyProviderLoad(tool.get(), verr)) {
      syst_param_id_offset += tool->GetSystMetaData().size();
      verified.push_back(std::move(tool));
    } else {
      err_out = "verify failed: " + verr;
      return std::vector<std::unique_ptr<IGENIESystProvider_tool>>{};
    }
  }
  return verified;
}

// ===== main ================================================================
int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  // Resolve fcl-dir default
  if (cliopts::fcl_dir.empty()) {
    char const *nusyst_env = std::getenv("NUSYST");
    if (nusyst_env) cliopts::fcl_dir = std::string(nusyst_env) + "/fcl";
  }

  // Collect all successfully-loaded providers here
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> all_providers;
  // Running paramId offset, threaded through every loader call so that
  // parameter IDs remain unique across providers in the merged output.
  paramId_t syst_param_id_offset = 0;

  // ----- GENIE RW mode -----
  if (cliopts::mode == "genierw" || cliopts::mode == "all") {
    std::cerr << "=== GENIE Reweight ===" << std::endl;
    std::vector<std::string> included, skipped;
    fhicl::ParameterSet grw_ps = BuildGENIEReWeightToolConfig(included, skipped);
    std::cerr << "Enumerated " << (included.size() + skipped.size())
              << " GSyst_t entries (" << included.size() << " active, "
              << skipped.size() << " skipped)" << std::endl;
    if (!skipped.empty()) {
      std::cerr << "Skipped GENIE RW dials:" << std::endl;
      for (auto &s : skipped) std::cerr << "    " << s << std::endl;
    }

    std::string err;
    auto tools = TryLoadInMemoryGENIERW(grw_ps, syst_param_id_offset, err);
    if (tools.empty()) {
      std::cerr << "[ERROR]: GENIE Reweight provider failed to instantiate:\n"
                << "  " << err << std::endl;
      std::cerr << "  Hint: ensure GENIE_XSEC_TUNE is set and try --include-zexp\n"
                << "        if the dipole CCQE FF dials are not the desired set."
                << std::endl;
    } else {
      for (auto &t : tools) {
        std::cerr << "  loaded: " << t->GetFullyQualifiedName() << " ("
                  << t->GetSystMetaData().size() << " params)" << std::endl;
        all_providers.push_back(std::move(t));
      }
    }
  }

  // ----- providers mode -----
  if (cliopts::mode == "providers" || cliopts::mode == "all") {
    std::cerr << "\n=== Standalone providers ===" << std::endl;
    if (cliopts::fcl_dir.empty()) {
      std::cerr << "[ERROR]: --fcl-dir not set and $NUSYST/fcl unavailable; "
                << "skipping providers." << std::endl;
    } else {
      std::cerr << "Scanning " << cliopts::fcl_dir << std::endl;
      std::vector<std::string> fcls = ListFcls(cliopts::fcl_dir);
      std::cerr << "Found " << fcls.size() << " .fcl files" << std::endl;

      for (auto &path : fcls) {
        std::string base = path.substr(path.find_last_of('/') + 1);

        // Skip user-specified files
        bool skip = false;
        for (auto &s : cliopts::skip_providers) {
          if (base == s) { skip = true; break; }
        }
        if (skip) {
          std::cerr << "  [SKIP] " << base << " (user --skip)" << std::endl;
          continue;
        }

        // Skip non-toolconfig files (e.g. paramHeader_FSI.fcl is a generated config)
        if (!IsToolConfigFile(path)) {
          std::cerr << "  [SKIP] " << base << " (not a tool config: no syst_providers key)" << std::endl;
          continue;
        }

        std::string err;
        auto tools = TryLoadToolConfig(path, syst_param_id_offset, err);
        if (tools.empty()) {
          std::cerr << "  [FAIL] " << base << "\n         " << err << std::endl;
          continue;
        }
        std::cerr << "  [OK]   " << base;
        for (auto &t : tools) {
          std::cerr << "  -> " << t->GetFullyQualifiedName() << " ("
                    << t->GetSystMetaData().size() << " params)";
          all_providers.push_back(std::move(t));
        }
        std::cerr << std::endl;
      }
    }
  }

  if (all_providers.empty()) {
    std::cerr << "\n[ERROR]: No providers loaded successfully. Nothing to write."
              << std::endl;
    return 2;
  }

  // ----- Build the merged generated config -----
  fhicl::ParameterSet out_ps;
  std::vector<std::string> providerNames;
  for (auto &prov : all_providers) {
    if (!systtools::Validate(prov->GetSystMetaData(), false)) {
      std::cerr << "[WARN]: Provider " << prov->GetFullyQualifiedName()
                << " failed Validate(); skipping." << std::endl;
      continue;
    }
    fhicl::ParameterSet tool_ps = prov->GetParameterHeadersDocument();
    out_ps.put(prov->GetFullyQualifiedName(), tool_ps);
    providerNames.push_back(prov->GetFullyQualifiedName());
  }
  out_ps.put("syst_providers", providerNames);

  fhicl::ParameterSet wrapped_out_ps;
  wrapped_out_ps.put("generated_systematic_provider_configuration", out_ps);

  std::ostream *os(nullptr);
  std::ofstream fs;
  if (cliopts::outputfile.size()) {
    fs.open(cliopts::outputfile);
    if (!fs.is_open()) {
      std::cerr << "[ERROR]: Failed to open " << cliopts::outputfile << std::endl;
      return 1;
    }
    os = &fs;
  } else {
    os = &std::cout;
  }
  (*os) << wrapped_out_ps.to_indented_string() << std::endl;
  if (cliopts::outputfile.size()) fs.close();

  std::cerr << "\n=== Summary ===" << std::endl;
  std::cerr << "Wrote " << providerNames.size() << " providers" << std::endl;
  if (cliopts::outputfile.size())
    std::cerr << "Output: " << cliopts::outputfile << std::endl;

  return 0;
}
