// DeclaredDialTestNuSyst.cxx
// Prints an inventory of all dials declared in a nusystematics fhicl config,
// with their metadata (CV value, one-sigma shifts, variations, flags).
// Inspired by NIWGReWeight's DeclaredDialTest.

#include "systematicstools/interface/SystParamHeader.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/response_helper.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace systtools;
using namespace nusyst;

namespace cliopts {
std::string fclname;
std::vector<std::string> parameters;
bool verbose = false;
}

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Required:\n"
    "    -c <config.fcl>      FHiCL config file\n\n"
    "  Optional:\n"
    "    -p <par1,par2,...>    Filter dials by name substring\n"
    "    --verbose             Print paramVariations, validity range, opts\n"
    "    -?|--help             Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "-c") cliopts::fclname = argv[++opt];
    else if (s == "--verbose") cliopts::verbose = true;
    else if (s == "-p") {
      std::string tok; std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::parameters.push_back(tok);
    } else {
      std::cout << "[ERROR]: Unknown option: " << s << std::endl;
      SayUsage(argv); exit(1);
    }
    opt++;
  }
}

bool ParamSelected(const std::string &name) {
  if (cliopts::parameters.empty()) return true;
  for (auto &sel : cliopts::parameters)
    if (name.find(sel) != std::string::npos) return true;
  return false;
}

std::string FmtDouble(double v, int width = 8) {
  if (v == kDefaultDouble) return std::string(width - 3, ' ') + "---";
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(3) << std::setw(width) << v;
  return ss.str();
}

std::string FmtVariations(const std::vector<double> &vars) {
  if (vars.empty()) return "(none)";
  std::ostringstream ss;
  ss << "[";
  for (size_t i = 0; i < vars.size(); ++i) {
    if (i > 0) ss << ", ";
    ss << std::fixed << std::setprecision(3) << vars[i];
  }
  ss << "]";
  return ss.str();
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  if (cliopts::fclname.empty()) {
    std::cout << "[ERROR]: -c <config.fcl> is required." << std::endl;
    SayUsage(argv);
    return 1;
  }

  response_helper phh(cliopts::fclname);

  // Collect rows first so we can align columns
  struct Row {
    paramId_t pid;
    std::string provider;
    std::string name;
    double cv;
    double sigma_up, sigma_down;
    int nvar;
    bool splineable, correction, natural, weight_syst;
    std::vector<double> variations;
    std::array<double, 2> validity;
    std::vector<std::string> opts;
  };
  std::vector<Row> rows;

  for (paramId_t pid : phh.GetParameters()) {
    SystParamHeader const &hdr = phh.GetHeader(pid);
    if (!ParamSelected(hdr.prettyName)) continue;

    // Find provider name (pattern from DumpConfiguredTweaksNuSyst)
    std::string provname = "?";
    for (size_t si = 0; si < phh.GetSystProvider().size(); ++si) {
      if (phh.GetSystProvider()[si]->ParamIsHandled(pid)) {
        provname = phh.GetSystProvider()[si]->GetFullyQualifiedName();
        break;
      }
    }

    Row r;
    r.pid = pid;
    r.provider = provname;
    r.name = hdr.prettyName;
    r.cv = hdr.centralParamValue;
    r.sigma_up = hdr.oneSigmaShifts[1];
    r.sigma_down = hdr.oneSigmaShifts[0];
    r.nvar = (int)hdr.paramVariations.size();
    r.splineable = hdr.isSplineable;
    r.correction = hdr.isCorrection;
    r.natural = hdr.unitsAreNatural;
    r.weight_syst = hdr.isWeightSystematicVariation;
    r.variations = hdr.paramVariations;
    r.validity = hdr.paramValidityRange;
    r.opts = hdr.opts;
    rows.push_back(std::move(r));
  }

  if (rows.empty()) {
    std::cout << "[INFO]: No dials matched the selection." << std::endl;
    return 0;
  }

  // Column widths
  size_t w_id = 4, w_prov = 10, w_name = 10;
  for (auto &r : rows) {
    w_id = std::max(w_id, std::to_string(r.pid).size());
    w_prov = std::max(w_prov, r.provider.size());
    w_name = std::max(w_name, r.name.size());
  }
  w_id += 2; w_prov += 2; w_name += 2;

  // Header
  std::cout << "\n";
  std::cout << "=== Dial inventory from " << cliopts::fclname << " ===\n";
  std::cout << "Total dials: " << rows.size() << "\n\n";

  std::cout << std::left
            << std::setw(w_id) << "ID"
            << std::setw(w_prov) << "Provider"
            << std::setw(w_name) << "Name"
            << std::setw(10) << "CV"
            << std::setw(10) << "+1sigma"
            << std::setw(10) << "-1sigma"
            << std::setw(6) << "Nvar"
            << std::setw(7) << "spline"
            << std::setw(7) << "corr"
            << std::setw(7) << "units"
            << "\n";

  size_t total_w = w_id + w_prov + w_name + 10 + 10 + 10 + 6 + 7 + 7 + 7;
  std::cout << std::string(total_w, '-') << "\n";

  for (auto &r : rows) {
    std::cout << std::left
              << std::setw(w_id) << r.pid
              << std::setw(w_prov) << r.provider
              << std::setw(w_name) << r.name
              << std::setw(10) << FmtDouble(r.cv)
              << std::setw(10) << FmtDouble(r.sigma_up)
              << std::setw(10) << FmtDouble(r.sigma_down)
              << std::setw(6) << r.nvar
              << std::setw(7) << (r.splineable ? "yes" : "no")
              << std::setw(7) << (r.correction ? "yes" : "no")
              << std::setw(7) << (r.natural ? "nat" : "sigma")
              << "\n";

    if (cliopts::verbose) {
      std::cout << "    paramVariations: " << FmtVariations(r.variations) << "\n";
      std::cout << "    validityRange:   ["
                << FmtDouble(r.validity[0], 6) << ", "
                << FmtDouble(r.validity[1], 6) << "]\n";
      std::cout << "    weightSyst:      " << (r.weight_syst ? "yes" : "no (property shift)") << "\n";
      if (!r.opts.empty()) {
        std::cout << "    opts:            ";
        for (size_t i = 0; i < r.opts.size(); ++i) {
          if (i > 0) std::cout << ", ";
          std::cout << r.opts[i];
        }
        std::cout << "\n";
      }
      std::cout << "\n";
    }
  }

  std::cout << std::endl;
  return 0;
}
