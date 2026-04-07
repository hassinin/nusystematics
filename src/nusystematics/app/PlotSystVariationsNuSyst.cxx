// PlotSystVariationsNuSyst.cxx
// Plots systematic variations for nusystematics parameters,
// separated by interaction channel (QE, MEC, RES, DIS, Other)
// and by CC/NC, neutrino species, and target.
// Produces differential cross section plots with ratio panels,
// a summary page, and a shared legend.

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"
#include "nusystematics/utility/response_helper.hh"

#include "fhiclcpp/ParameterSet.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Conventions/Units.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObjString.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace systtools;
using namespace nusyst;

// ===== Predefined variables ================================================
struct VarDef {
  std::string branch;
  std::string label;      // x-axis label
  std::string diffLabel;  // y-axis for diff xsec
  double xmin, xmax;
};

static const std::map<std::string, VarDef> kPredefinedVars = {
  {"Enu",      {"Enu_true",          "E_{#nu} (GeV)",      "d#sigma/dE_{#nu}",          0, 10}},
  {"q0",       {"q0",                "q_{0} (GeV)",        "d#sigma/dq_{0}",            0, 5}},
  {"Q2",       {"Q2",                "Q^{2} (GeV^{2})",    "d#sigma/dQ^{2}",            0, 5}},
  {"q3",       {"q3",                "|#vec{q}| (GeV)",    "d#sigma/d|#vec{q}|",        0, 5}},
  {"W",        {"w",                 "W (GeV)",            "d#sigma/dW",                0, 3}},
  {"plep",     {"plep",              "p_{lep} (GeV)",      "d#sigma/dp_{lep}",          0, 5}},
  {"coslep",   {"coslep",            "cos#theta_{lep}",    "d#sigma/dcos#theta_{lep}",  -1, 1}},
  {"EAvail",   {"EAvail_GeV",        "E_{avail} (GeV)",    "d#sigma/dE_{avail}",        0, 3}},
  {"Tp",       {"leading_proton_KE", "T_{p} (GeV)",        "d#sigma/dT_{p}",            0, 1.5}},
  {"Emiss",    {"Emiss",             "E_{miss} (GeV)",     "d#sigma/dE_{miss}",         0, 0.5}},
  {"pmiss",    {"pmiss",             "p_{miss} (GeV)",     "d#sigma/dp_{miss}",         0, 1}},
  {"Ereco",    {"Ereco_cal",         "E_{reco}^{cal} (GeV)","d#sigma/dE_{reco}",        0, 10}},
  {"nproton",  {"nproton",           "N_{proton}",         "d#sigma/dN_{p}",            0, 6}},
  {"npip",     {"npip",              "N_{#pi^{+}}",        "d#sigma/dN_{#pi^{+}}",      0, 5}},
  {"npim",     {"npim",              "N_{#pi^{-}}",        "d#sigma/dN_{#pi^{-}}",      0, 5}},
  {"npi0",     {"npi0",              "N_{#pi^{0}}",        "d#sigma/dN_{#pi^{0}}",      0, 5}},
  {"nneutron", {"nneutron",          "N_{neutron}",        "d#sigma/dN_{n}",            0, 6}},
};

VarDef ResolveVar(const std::string &name) {
  auto it = kPredefinedVars.find(name);
  if (it != kPredefinedVars.end()) return it->second;
  return {name, name, "d#sigma/d" + name, 0, 0};
}

// ===== Channel classification ==============================================
struct EventVars {
  double Enu, q0, Q2, q3, w, plep, coslep, EAvail, leading_proton_KE, xsec;
  double Emiss, pmiss, Ereco_cal;
  int nproton, npip, npim, npi0, nneutron;
  bool is_cc, is_qe, is_mec, is_res, is_dis;
  int nu_pdg, tgt_A, tgt_Z;
  bool has_coslep, has_leading_proton;
};

std::string NuName(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "numu" : "numubar";
    case 12: return pdg > 0 ? "nue"  : "nuebar";
    case 16: return pdg > 0 ? "nutau" : "nutaubar";
    default: return Form("nu%d", pdg);
  }
}

std::string NuLabel(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "#nu_{#mu}" : "#bar{#nu}_{#mu}";
    case 12: return pdg > 0 ? "#nu_{e}"   : "#bar{#nu}_{e}";
    case 16: return pdg > 0 ? "#nu_{#tau}": "#bar{#nu}_{#tau}";
    default: return Form("#nu_{%d}", pdg);
  }
}

std::string TargetName(int A) {
  switch (A) {
    case 40: return "Ar40";
    case 12: return "C12";
    case 16: return "O16";
    case 56: return "Fe56";
    case 1:  return "H";
    default: return Form("A%d", A);
  }
}

std::string TargetLabel(int A) {
  switch (A) {
    case 40: return "^{40}Ar";
    case 12: return "^{12}C";
    case 16: return "^{16}O";
    case 56: return "^{56}Fe";
    case 1:  return "H";
    default: return Form("A=%d", A);
  }
}

std::string InteractionType(const EventVars &ev) {
  if (ev.is_qe)  return "QE";
  if (ev.is_mec) return "MEC";
  if (ev.is_res) return "RES";
  if (ev.is_dis) return "DIS";
  return "Other";
}

std::string InteractionLabel(const std::string &type) {
  if (type == "MEC") return "2p2h";
  return type;
}

// Sort key: CC before NC, numu before nue, QE/MEC/RES/DIS/Other order
std::string ChannelSortKey(const std::string &channel) {
  std::string key;
  key += (channel.find("CC_") == 0) ? "0" : "1";
  if (channel.find("numu_")    != std::string::npos) key += "0";
  else if (channel.find("nue_") != std::string::npos) key += "1";
  else if (channel.find("nutau_") != std::string::npos) key += "2";
  else key += "3";
  if (channel.find("bar_") != std::string::npos) key += "1"; else key += "0";
  if (channel.find("_QE") != std::string::npos) key += "0";
  else if (channel.find("_MEC") != std::string::npos) key += "1";
  else if (channel.find("_RES") != std::string::npos) key += "2";
  else if (channel.find("_DIS") != std::string::npos) key += "3";
  else key += "4";
  return key + "_" + channel;
}

std::string MakeChannelKey(const EventVars &ev) {
  std::string cc = ev.is_cc ? "CC" : "NC";
  return cc + "_" + NuName(ev.nu_pdg) + "_" + TargetName(ev.tgt_A) + "_" + InteractionType(ev);
}

std::string MakeChannelLabel(const std::string &key) {
  if (key == "Total") return "Total (all channels)";
  if (key == "Total_CC") return "Total CC";
  if (key == "Total_NC") return "Total NC";
  // Parse key: CC_numu_Ar40_QE
  std::istringstream ss(key);
  std::string cc, nu, tgt, inttype;
  std::getline(ss, cc, '_');
  // nu name may contain "bar" so read until target
  std::string rest;
  std::getline(ss, rest);
  // Find last two underscore-separated tokens as target and interaction
  size_t last_us = rest.rfind('_');
  inttype = rest.substr(last_us + 1);
  rest = rest.substr(0, last_us);
  size_t second_last = rest.rfind('_');
  tgt = rest.substr(second_last + 1);
  std::string nuname = rest.substr(0, second_last);

  // Build label
  int pdg = 0;
  if (nuname == "numu") pdg = 14;
  else if (nuname == "numubar") pdg = -14;
  else if (nuname == "nue") pdg = 12;
  else if (nuname == "nuebar") pdg = -12;
  else if (nuname == "nutau") pdg = 16;
  else if (nuname == "nutaubar") pdg = -16;

  int A = 0;
  if (tgt == "Ar40") A = 40;
  else if (tgt == "C12") A = 12;
  else if (tgt == "O16") A = 16;
  else if (tgt == "Fe56") A = 56;
  else if (tgt == "H") A = 1;
  else sscanf(tgt.c_str(), "A%d", &A);

  return cc + " " + (pdg ? NuLabel(pdg) : nuname) + " " +
         (A ? TargetLabel(A) : tgt) + " " + InteractionLabel(inttype);
}

// ===== Histogram bookkeeping ===============================================
struct SystHists {
  TH1D *h_cv = nullptr;
  std::vector<TH1D *> h_var;
  std::vector<double> tweakvals;
};

struct ParamMeta {
  std::string fullname;
  int ntweaks;
  std::vector<double> tweakvalues;
};

// ===== CLI =================================================================
namespace cliopts {
std::string input_file;
std::string fclname;
std::string output_base = "syst_variations";
std::string treename = "events";
std::string genie_branch_name = "gmcrec";
std::vector<std::string> variables;
std::vector<std::string> parameters;
size_t NMax = std::numeric_limits<size_t>::max();
int nbins = 50;
bool make_pdf = true;
bool make_root = true;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Required:\n"
    "    -i <file.root>       Input ROOT file (ghep or flat tree)\n\n"
    "  Optional:\n"
    "    -c <config.fcl>      FHiCL config (required for GHEP mode)\n"
    "    -o <output>          Output basename (default: syst_variations)\n"
    "    -v <var1,var2,...>    Variables to plot (default: all predefined)\n"
    "                         Predefined: Enu,q0,Q2,q3,W,plep,coslep,EAvail,Tp\n"
    "    -p <par1,par2,...>    Systematic parameters to plot (default: all)\n"
    "    -t <treename>        Tree name (default: events)\n"
    "    -b <branch>          NtpMCEventRecord branch (default: gmcrec)\n"
    "    -N <NMax>            Max events to process\n"
    "    --nbins <n>          Number of bins (default: 50)\n"
    "    --no-pdf / --no-root Skip output\n"
    "    -?|--help            Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "-i") cliopts::input_file = argv[++opt];
    else if (s == "-c") cliopts::fclname = argv[++opt];
    else if (s == "-o") cliopts::output_base = argv[++opt];
    else if (s == "-t") cliopts::treename = argv[++opt];
    else if (s == "-b") cliopts::genie_branch_name = argv[++opt];
    else if (s == "-N") cliopts::NMax = str2T<size_t>(argv[++opt]);
    else if (s == "--nbins") cliopts::nbins = str2T<int>(argv[++opt]);
    else if (s == "--no-pdf") cliopts::make_pdf = false;
    else if (s == "--no-root") cliopts::make_root = false;
    else if (s == "-v") {
      std::string tok; std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::variables.push_back(tok);
    } else if (s == "-p") {
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

// ===== Event counting for summary ==========================================
struct EventCounts {
  size_t total = 0;
  std::map<int, size_t> by_nu_pdg;
  std::map<int, size_t> by_tgt_A;
  std::map<std::string, size_t> by_channel; // full channel key
  std::map<std::string, size_t> by_current; // "CC" or "NC"
  std::map<std::string, size_t> by_inttype; // "CC_QE", "NC_RES", etc.

  void Count(const EventVars &ev) {
    total++;
    by_nu_pdg[ev.nu_pdg]++;
    by_tgt_A[ev.tgt_A]++;
    std::string cc = ev.is_cc ? "CC" : "NC";
    by_current[cc]++;
    by_inttype[cc + "_" + InteractionType(ev)]++;
    by_channel[MakeChannelKey(ev)]++;
  }
};

// ===== Helpers =============================================================
std::pair<double, double> AutoRange(TTree *tree, const std::string &branch) {
  double val;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(branch.c_str(), 1);
  if (!tree->GetBranch(branch.c_str())) { tree->SetBranchStatus("*", 1); return {0, 1}; }
  tree->SetBranchAddress(branch.c_str(), &val);
  double vmin = 1e30, vmax = -1e30;
  Long64_t n = std::min((Long64_t)10000, tree->GetEntries());
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    if (std::isfinite(val)) { vmin = std::min(vmin, val); vmax = std::max(vmax, val); }
  }
  tree->SetBranchStatus("*", 1);
  tree->ResetBranchAddresses();
  if (vmin >= vmax) return {0, 1};
  double pad = 0.05 * (vmax - vmin);
  return {vmin - pad, vmax + pad};
}

void BookHistograms(const std::string &channel, const ParamMeta &pm,
                    const std::vector<VarDef> &vars,
                    std::map<std::string, SystHists> &var_hists,
                    TTree *tree_for_autorange = nullptr) {
  for (auto &vd : vars) {
    double xmin = vd.xmin, xmax = vd.xmax;
    if (xmin == 0 && xmax == 0 && tree_for_autorange) {
      auto r = AutoRange(tree_for_autorange, vd.branch);
      xmin = r.first; xmax = r.second;
    }
    if (xmin == 0 && xmax == 0) { xmin = 0; xmax = 5; }

    std::string tag = channel + "_" + pm.fullname + "_" + vd.branch;
    SystHists sh;
    sh.h_cv = new TH1D(("h_cv_" + tag).c_str(), "", cliopts::nbins, xmin, xmax);
    sh.h_cv->Sumw2();
    sh.tweakvals = pm.tweakvalues;
    for (int i = 0; i < pm.ntweaks; ++i) {
      TH1D *hv = new TH1D(Form("h_var%d_%s", i, tag.c_str()), "", cliopts::nbins, xmin, xmax);
      hv->Sumw2();
      sh.h_var.push_back(hv);
    }
    var_hists[vd.branch] = sh;
  }
}

double GetVarValue(const EventVars &ev, const std::string &branch) {
  if (branch == "Enu_true") return ev.Enu;
  if (branch == "q0") return ev.q0;
  if (branch == "Q2") return ev.Q2;
  if (branch == "q3") return ev.q3;
  if (branch == "w") return ev.w;
  if (branch == "plep") return ev.plep;
  if (branch == "coslep") return ev.has_coslep ? ev.coslep : -999;
  if (branch == "EAvail_GeV") return ev.EAvail;
  if (branch == "leading_proton_KE") return ev.has_leading_proton ? ev.leading_proton_KE : -999;
  if (branch == "Emiss") return ev.Emiss;
  if (branch == "pmiss") return ev.pmiss;
  if (branch == "Ereco_cal") return ev.Ereco_cal;
  if (branch == "nproton") return (double)ev.nproton;
  if (branch == "npip") return (double)ev.npip;
  if (branch == "npim") return (double)ev.npim;
  if (branch == "npi0") return (double)ev.npi0;
  if (branch == "nneutron") return (double)ev.nneutron;
  return -999;
}

// ===== Histogram map type ==================================================
// allhists[channel][param][var] -> SystHists
using HistMap = std::map<std::string, std::map<std::string, std::map<std::string, SystHists>>>;

// ===== FLAT TREE MODE ======================================================
void FillFromFlatTree(
    TTree *tree, TTree *meta,
    const std::vector<VarDef> &vars,
    const std::vector<ParamMeta> &params,
    HistMap &allhists,
    EventCounts &counts) {

  // Branch addresses for event classification
  double xsec = 0;
  int nu_pdg = 0, tgt_A = 0, tgt_Z = 0;
  bool is_cc = false, is_qe = false, is_mec = false, is_res = false, is_dis = false;
  double q0 = 0, Q2 = 0, q3 = 0, w = 0, plep = 0, Enu_true = 0, EAvail_GeV = 0;
  double Emiss = 0, pmiss = 0;
  std::vector<double> *fsprotons_KE = nullptr;
  std::vector<int> *fsi_pdgs = nullptr;

  tree->SetBranchAddress("xsec", &xsec);
  tree->SetBranchAddress("nu_pdg", &nu_pdg);
  tree->SetBranchAddress("tgt_A", &tgt_A);
  tree->SetBranchAddress("tgt_Z", &tgt_Z);
  tree->SetBranchAddress("is_cc", &is_cc);
  tree->SetBranchAddress("is_qe", &is_qe);
  tree->SetBranchAddress("is_mec", &is_mec);
  tree->SetBranchAddress("is_res", &is_res);
  tree->SetBranchAddress("is_dis", &is_dis);
  tree->SetBranchAddress("q0", &q0);
  tree->SetBranchAddress("Q2", &Q2);
  tree->SetBranchAddress("q3", &q3);
  tree->SetBranchAddress("w", &w);
  tree->SetBranchAddress("plep", &plep);
  tree->SetBranchAddress("Enu_true", &Enu_true);
  tree->SetBranchAddress("EAvail_GeV", &EAvail_GeV);

  bool has_Emiss = tree->GetBranch("Emiss") != nullptr;
  bool has_pmiss = tree->GetBranch("pmiss") != nullptr;
  if (has_Emiss) tree->SetBranchAddress("Emiss", &Emiss);
  if (has_pmiss) tree->SetBranchAddress("pmiss", &pmiss);

  bool has_proton_branch = tree->GetBranch("fsprotons_KE") != nullptr;
  if (has_proton_branch) tree->SetBranchAddress("fsprotons_KE", &fsprotons_KE);

  bool has_fsi_pdgs = tree->GetBranch("fsi_pdgs") != nullptr;
  if (has_fsi_pdgs) tree->SetBranchAddress("fsi_pdgs", &fsi_pdgs);

  // Tweak response branches
  struct PBranch {
    std::string fullname;
    int ntweaks;
    std::vector<double> responses;
    double cv_weight;
  };
  std::vector<PBranch> pbranches;
  for (auto &pm : params) {
    PBranch pb;
    pb.fullname = pm.fullname;
    pb.ntweaks = pm.ntweaks;
    pb.responses.resize(pm.ntweaks, 1.0);
    pb.cv_weight = 1.0;
    if (tree->GetBranch(("tweak_responses_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("tweak_responses_" + pm.fullname).c_str(), pb.responses.data());
    if (tree->GetBranch(("paramCVWeight_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("paramCVWeight_" + pm.fullname).c_str(), &pb.cv_weight);
    pbranches.push_back(std::move(pb));
  }

  Long64_t nentries = std::min((Long64_t)cliopts::NMax, tree->GetEntries());
  Long64_t nshout = std::max(nentries / 20, (Long64_t)1);

  for (Long64_t ev = 0; ev < nentries; ++ev) {
    tree->GetEntry(ev);
    if (!(ev % nshout))
      std::cout << "\rProcessing event " << ev << "/" << nentries << std::flush;

    EventVars evars;
    evars.Enu = Enu_true; evars.q0 = q0; evars.Q2 = Q2; evars.q3 = q3;
    evars.w = w; evars.plep = plep; evars.EAvail = EAvail_GeV; evars.xsec = xsec;
    evars.is_cc = is_cc; evars.is_qe = is_qe; evars.is_mec = is_mec;
    evars.is_res = is_res; evars.is_dis = is_dis;
    evars.nu_pdg = nu_pdg; evars.tgt_A = tgt_A; evars.tgt_Z = tgt_Z;
    evars.has_coslep = false; evars.coslep = -999;
    evars.has_leading_proton = has_proton_branch && fsprotons_KE && !fsprotons_KE->empty();
    evars.leading_proton_KE = evars.has_leading_proton ? (*fsprotons_KE)[0] : -999;
    evars.Emiss = has_Emiss ? Emiss : -999;
    evars.pmiss = has_pmiss ? pmiss : -999;

    // Calorimetric reco energy: Eavail + lepton energy (simple calorimetric)
    double Elep = std::sqrt(plep * plep + 0.10566 * 0.10566); // assume muon mass
    evars.Ereco_cal = EAvail_GeV + Elep;

    // Particle multiplicities from flat tree
    evars.nproton = has_proton_branch && fsprotons_KE ? (int)fsprotons_KE->size() : 0;
    evars.npip = 0; evars.npim = 0; evars.npi0 = 0; evars.nneutron = 0;
    // fsi_pdgs is the list of FS hadron PDGs — but the flat tree stores these
    // as intranuclear particles, not final state. We don't have FS particle
    // multiplicities directly for pions/neutrons in the flat tree.
    // For now set to 0; GHEP mode will compute them properly.

    counts.Count(evars);
    std::string chkey = MakeChannelKey(evars);

    for (size_t ip = 0; ip < pbranches.size(); ++ip) {
      auto &pb = pbranches[ip];
      auto &pm = params[ip];

      // Book histograms on first encounter of this channel
      if (allhists[chkey].find(pm.fullname) == allhists[chkey].end()) {
        BookHistograms(chkey, pm, vars, allhists[chkey][pm.fullname], tree);
      }

      auto &var_hists = allhists[chkey][pm.fullname];
      for (auto &vd : vars) {
        double x = GetVarValue(evars, vd.branch);
        if (x <= -999) continue;
        auto vit = var_hists.find(vd.branch);
        if (vit == var_hists.end()) continue;
        auto &sh = vit->second;

        sh.h_cv->Fill(x, xsec * pb.cv_weight);
        for (int i = 0; i < pb.ntweaks && i < (int)sh.h_var.size(); ++i)
          sh.h_var[i]->Fill(x, xsec * pb.cv_weight * pb.responses[i]);
      }
    }
  }
  std::cout << "\rProcessed " << nentries << " events.          " << std::endl;
  tree->ResetBranchAddresses();
}

// ===== GHEP MODE ===========================================================
void FillFromGHEP(
    const std::string &input, const std::string &fclname,
    const std::vector<VarDef> &vars,
    HistMap &allhists,
    EventCounts &counts) {

  response_helper phh(fclname);

  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(input.c_str())) {
    std::cout << "[ERROR]: No gtree found in " << input << std::endl;
    return;
  }
  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  gevs->SetBranchAddress(cliopts::genie_branch_name.c_str(), &GenieNtpl);

  // Build param metadata
  struct GPInfo { paramId_t pid; ParamMeta meta; };
  std::vector<GPInfo> gparams;
  for (paramId_t pid : phh.GetParameters()) {
    SystParamHeader const &hdr = phh.GetHeader(pid);
    if (hdr.isResponselessParam) continue;
    std::string provname;
    for (size_t si = 0; si < phh.GetSystProvider().size(); ++si)
      if (phh.GetSystProvider()[si]->ParamIsHandled(pid))
        { provname = phh.GetSystProvider()[si]->GetFullyQualifiedName(); break; }
    std::string fullname = provname + "_" + hdr.prettyName;
    if (!ParamSelected(fullname)) continue;
    GPInfo gi;
    gi.pid = pid;
    gi.meta.fullname = fullname;
    gi.meta.ntweaks = hdr.isCorrection ? 1 : (int)hdr.paramVariations.size();
    gi.meta.tweakvalues = hdr.paramVariations;
    gparams.push_back(gi);
  }
  if (gparams.empty()) { std::cout << "[WARN]: No matching parameters." << std::endl; delete gevs; return; }

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");

  Long64_t nentries = std::min((Long64_t)cliopts::NMax, gevs->GetEntries());
  Long64_t nshout = std::max(nentries / 20, (Long64_t)1);

  for (Long64_t ev_it = 0; ev_it < nentries; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;
    if (!(ev_it % nshout))
      std::cout << "\rGHEP event " << ev_it << "/" << nentries << std::flush;

    genie::GHepParticle *FSLep = GenieGHep.FinalStatePrimaryLepton();
    genie::GHepParticle *ISLep = GenieGHep.Probe();
    if (!FSLep || !ISLep) { GenieNtpl->Clear(); continue; }

    TLorentzVector FSLepP4 = *FSLep->P4();
    TLorentzVector ISLepP4 = *ISLep->P4();
    TLorentzVector emTransfer = ISLepP4 - FSLepP4;

    EventVars evars;
    evars.xsec = GenieGHep.XSec() / genie::units::cm2;
    evars.Enu = ISLepP4.E();
    evars.q0 = emTransfer.E();
    evars.Q2 = -emTransfer.Mag2();
    evars.q3 = emTransfer.Vect().Mag();
    double MN = 0.93827203;
    double W2 = MN*MN + 2.0*MN*evars.q0 - evars.Q2;
    evars.w = (W2 > 0) ? std::sqrt(W2) : -1;
    evars.plep = FSLepP4.Vect().Mag();
    evars.has_coslep = true;
    evars.coslep = FSLepP4.CosTheta();
    evars.EAvail = GetErecoil_MINERvA_LowRecoil(GenieGHep);

    // FS particle loop: proton KE, multiplicities
    evars.has_leading_proton = false;
    evars.leading_proton_KE = -999;
    evars.nproton = 0; evars.npip = 0; evars.npim = 0; evars.npi0 = 0; evars.nneutron = 0;
    genie::GHepParticle *p = nullptr;
    TIter iter(&GenieGHep);
    double max_pKE = -1;
    while ((p = dynamic_cast<genie::GHepParticle*>(iter.Next()))) {
      if (p->Status() != genie::kIStStableFinalState) continue;
      int pdgc = p->Pdg();
      if (pdgc == 2212) {
        evars.nproton++;
        double ke = p->KinE();
        if (ke > max_pKE) { max_pKE = ke; evars.has_leading_proton = true; evars.leading_proton_KE = ke; }
      }
      else if (pdgc ==  211) evars.npip++;
      else if (pdgc == -211) evars.npim++;
      else if (pdgc ==  111) evars.npi0++;
      else if (pdgc == 2112) evars.nneutron++;
    }

    evars.Emiss = GetEmiss(GenieGHep, false);
    evars.pmiss = GetPmiss(GenieGHep, false);
    double Elep = std::sqrt(evars.plep * evars.plep + 0.10566 * 0.10566);
    evars.Ereco_cal = evars.EAvail + Elep;

    evars.is_cc = GenieGHep.Summary()->ProcInfo().IsWeakCC();
    evars.is_qe = GenieGHep.Summary()->ProcInfo().IsQuasiElastic();
    evars.is_mec = GenieGHep.Summary()->ProcInfo().IsMEC();
    evars.is_res = GenieGHep.Summary()->ProcInfo().IsResonant();
    evars.is_dis = GenieGHep.Summary()->ProcInfo().IsDeepInelastic();
    evars.nu_pdg = ISLep->Pdg();
    evars.tgt_A = GenieGHep.TargetNucleus() ? GenieGHep.TargetNucleus()->A() : 1;
    evars.tgt_Z = GenieGHep.TargetNucleus() ? GenieGHep.TargetNucleus()->Z() : 1;

    counts.Count(evars);
    std::string chkey = MakeChannelKey(evars);

    event_unit_response_w_cv_t resp = phh.GetEventVariationAndCVResponse(GenieGHep);

    for (auto &gp : gparams) {
      if (allhists[chkey].find(gp.meta.fullname) == allhists[chkey].end())
        BookHistograms(chkey, gp.meta, vars, allhists[chkey][gp.meta.fullname]);

      size_t ridx = GetParamContainerIndex(resp, gp.pid);
      double cv_w = 1.0;
      std::vector<double> responses;
      if (ridx != kParamUnhandled<size_t>) {
        cv_w = resp[ridx].CV_response;
        responses = resp[ridx].responses;
      }

      auto &var_hists = allhists[chkey][gp.meta.fullname];
      for (auto &vd : vars) {
        double x = GetVarValue(evars, vd.branch);
        if (x <= -999) continue;
        auto vit = var_hists.find(vd.branch);
        if (vit == var_hists.end()) continue;
        auto &sh = vit->second;
        sh.h_cv->Fill(x, evars.xsec * cv_w);
        for (int i = 0; i < (int)responses.size() && i < (int)sh.h_var.size(); ++i)
          sh.h_var[i]->Fill(x, evars.xsec * cv_w * responses[i]);
      }
    }
    GenieNtpl->Clear();
  }
  std::cout << "\rProcessed " << nentries << " GHEP events.          " << std::endl;
  delete gevs;
}

// ===== Scale histograms to differential xsec ===============================
void ScaleToDiffXsec(HistMap &allhists) {
  for (auto &[ch, param_map] : allhists)
    for (auto &[par, var_map] : param_map)
      for (auto &[var, sh] : var_map) {
        sh.h_cv->Scale(1.0, "width");
        for (auto *h : sh.h_var) h->Scale(1.0, "width");
      }
}

// ===== Remove empty channels ===============================================
void PruneEmpty(HistMap &allhists) {
  std::vector<std::string> empty_channels;
  for (auto &[ch, param_map] : allhists) {
    bool has_entries = false;
    for (auto &[par, var_map] : param_map)
      for (auto &[var, sh] : var_map)
        if (sh.h_cv->GetEntries() > 0) has_entries = true;
    if (!has_entries) empty_channels.push_back(ch);
  }
  for (auto &ch : empty_channels) allhists.erase(ch);
}

// ===== Color palette =======================================================
void GetVarColor(int ivar, int &color, int &style) {
  static const int palette[] = {
    kBlue+2, kBlue, kAzure+7, kGreen+3, kSpring+4,
    kOrange+7, kRed, kRed+2, kMagenta+2, kViolet+1
  };
  static const int npal = sizeof(palette) / sizeof(palette[0]);
  color = palette[ivar % npal];
  style = 1 + (ivar / npal);
}

// ===== Summary page ========================================================
void MakeSummaryPage(TCanvas *c, const EventCounts &counts) {
  c->Clear();
  c->cd();
  TPad *pad = new TPad("summary", "", 0, 0, 1, 1);
  pad->SetFillColor(kWhite);
  pad->Draw();
  pad->cd();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);

  double y = 0.92;
  double dy = 0.035;

  tex.SetTextSize(0.05);
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Event Summary"); y -= 1.5 * dy;
  tex.SetTextFont(42);
  tex.SetTextSize(0.032);

  tex.DrawLatex(0.05, y, Form("Total events processed: %zu", counts.total)); y -= 1.2 * dy;

  // Neutrino species
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Neutrino species:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[pdg, n] : counts.by_nu_pdg) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", NuName(pdg).c_str(), n));
    y -= dy;
  }

  // Targets
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Targets:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[A, n] : counts.by_tgt_A) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", TargetName(A).c_str(), n));
    y -= dy;
  }

  // By current
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "By current:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[cur, n] : counts.by_current) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", cur.c_str(), n));
    y -= dy;
  }

  // Interaction types
  double y_right = 0.92 - 1.5 * dy; // start on right column
  tex.SetTextFont(62);
  tex.DrawLatex(0.50, y_right, "Interaction breakdown:"); y_right -= dy;
  tex.SetTextFont(42);

  // Sort by key for consistent ordering
  std::vector<std::pair<std::string, size_t>> sorted_int(
      counts.by_inttype.begin(), counts.by_inttype.end());
  std::sort(sorted_int.begin(), sorted_int.end());
  for (auto &[key, n] : sorted_int) {
    tex.DrawLatex(0.53, y_right, Form("%-16s %zu", key.c_str(), n));
    y_right -= dy;
  }
}

// ===== PDF plotting (appends pages to already-open PDF) ====================
void MakePlots(const HistMap &allhists, const std::vector<VarDef> &vars) {
  std::string pdfname = cliopts::output_base + ".pdf";
  TCanvas *c = new TCanvas("cplots", "", 1400, 1000);

  // Sort channels: Total first, then Total_CC, Total_NC, then per-channel
  std::vector<std::string> channels;
  for (auto &[ch, _] : allhists) channels.push_back(ch);
  std::sort(channels.begin(), channels.end(),
            [](const std::string &a, const std::string &b) {
              // Priority: Total=0, Total_CC=1, Total_NC=2, rest sorted
              auto pri = [](const std::string &s) -> std::string {
                if (s == "Total") return "000";
                if (s == "Total_CC") return "001";
                if (s == "Total_NC") return "002";
                return "1" + ChannelSortKey(s);
              };
              return pri(a) < pri(b);
            });

  // Collect parameter names (ordered)
  std::set<std::string> param_set;
  for (auto &[ch, pm] : allhists)
    for (auto &[pname, _] : pm) param_set.insert(pname);
  std::vector<std::string> param_names(param_set.begin(), param_set.end());

  // One page per (parameter, channel)
  for (auto &parname : param_names) {
    for (auto &chkey : channels) {
      auto ch_it = allhists.find(chkey);
      if (ch_it == allhists.end()) continue;
      auto par_it = ch_it->second.find(parname);
      if (par_it == ch_it->second.end()) continue;
      auto &var_map = par_it->second;

      // Count valid (non-empty) variables for this page
      std::vector<const VarDef *> page_vars;
      for (auto &vd : vars) {
        auto vit = var_map.find(vd.branch);
        if (vit != var_map.end() && vit->second.h_cv->GetEntries() > 0)
          page_vars.push_back(&vd);
      }
      if (page_vars.empty()) continue;

      // Split into pages of max 6 plots each
      const int max_per_page = 6;
      int total_vars = (int)page_vars.size();
      int npages = (total_vars + max_per_page - 1) / max_per_page;

      for (int ipage = 0; ipage < npages; ++ipage) {
      int var_start = ipage * max_per_page;
      int var_end = std::min(var_start + max_per_page, total_vars);
      int nvars_this_page = var_end - var_start;

      int ncols = std::min(nvars_this_page, 3);
      int nrows = (nvars_this_page + ncols - 1) / ncols;

      c->Clear();
      c->SetFillColor(kWhite);

      // Title pad
      c->cd();
      TPad *titlepad = new TPad(Form("tp_%d", ipage), "", 0, 0.93, 1, 1.0);
      titlepad->SetFillColor(kWhite);
      titlepad->Draw(); titlepad->cd();
      TLatex title;
      title.SetNDC(); title.SetTextSize(0.42); title.SetTextAlign(22); title.SetTextFont(62);
      std::string page_title = parname + "  |  " + MakeChannelLabel(chkey);
      if (npages > 1) page_title += Form("  (%d/%d)", ipage + 1, npages);
      title.DrawLatex(0.5, 0.45, page_title.c_str());

      // Legend pad (horizontal, below title)
      c->cd();
      static int legpad_id = 0;
      TPad *legpad = new TPad(Form("lp_%d", legpad_id++), "", 0, 0.88, 1, 0.93);
      legpad->SetFillColor(kWhite);
      legpad->Draw(); legpad->cd();

      // Get variation values from the first variable's SystHists
      const SystHists &sh_first = var_map.begin()->second;
      int nvar = (int)sh_first.h_var.size();
      int nleg = nvar + 1;

      TLegend *leg = new TLegend(0.01, 0.05, 0.99, 0.95);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->SetTextSize(0.55); leg->SetNColumns(std::min(nleg, 10));

      // Dummy histograms for legend (unique names)
      static int hleg_uid = 0;
      TH1D *hleg_cv = new TH1D(Form("hleg_cv_%d", hleg_uid), "", 1, 0, 1);
      hleg_cv->SetLineColor(kBlack); hleg_cv->SetLineWidth(2);
      leg->AddEntry(hleg_cv, "CV", "l");

      for (int i = 0; i < nvar; ++i) {
        TH1D *hl = new TH1D(Form("hleg_%d_%d", hleg_uid, i), "", 1, 0, 1);
        int col_i, sty;
        GetVarColor(i, col_i, sty);
        hl->SetLineColor(col_i); hl->SetLineWidth(2); hl->SetLineStyle(sty);
        std::string lbl = (i < (int)sh_first.tweakvals.size()) ? Form("%.3g", sh_first.tweakvals[i]) : Form("var%d", i);
        leg->AddEntry(hl, lbl.c_str(), "l");
      }
      hleg_uid++;
      leg->Draw();

      // Grid of panels
      c->cd();
      double grid_top = 0.88, grid_bot = 0.0;
      double row_h = (grid_top - grid_bot) / nrows;
      double col_w = 1.0 / ncols;
      double ratio_frac = 0.30;

      for (int iv_abs = var_start; iv_abs < var_end; ++iv_abs) {
        int iv = iv_abs - var_start; // local index within this page
        const VarDef &vd = *page_vars[iv_abs];
        auto vit = var_map.find(vd.branch);
        if (vit == var_map.end()) continue;
        const SystHists &sh = vit->second;

        int col = iv % ncols;
        int row = iv / ncols;
        double x1 = col * col_w, x2 = (col + 1) * col_w;
        double y_top = grid_top - row * row_h;
        double y_bot = y_top - row_h;
        double y_split = y_bot + row_h * ratio_frac;

        // Main pad
        c->cd();
        TPad *pmain = new TPad(Form("m%d", iv), "", x1, y_split, x2, y_top);
        pmain->SetBottomMargin(0.01); pmain->SetTopMargin(0.04);
        pmain->SetLeftMargin(0.16); pmain->SetRightMargin(0.03);
        pmain->Draw(); pmain->cd();

        double ymax = sh.h_cv->GetMaximum();
        for (auto *h : sh.h_var) ymax = std::max(ymax, h->GetMaximum());

        TH1D *frame = (TH1D *)sh.h_cv->Clone(Form("f%d", iv));
        frame->Reset();
        frame->SetMaximum(ymax * 1.35); frame->SetMinimum(0);
        frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTickLength(0.03);
        frame->GetYaxis()->SetTitle(vd.diffLabel.c_str());
        frame->GetYaxis()->SetTitleSize(0.07); frame->GetYaxis()->SetTitleOffset(1.1);
        frame->GetYaxis()->SetLabelSize(0.06);
        frame->Draw("AXIS");

        sh.h_cv->SetLineColor(kBlack); sh.h_cv->SetLineWidth(2);
        sh.h_cv->Draw("HIST SAME");

        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          int col_i, sty;
          GetVarColor(i, col_i, sty);
          sh.h_var[i]->SetLineColor(col_i);
          sh.h_var[i]->SetLineWidth(2);
          sh.h_var[i]->SetLineStyle(sty);
          sh.h_var[i]->Draw("HIST SAME");
        }

        // Compute chi2/NDF and max fractional shift across all variations.
        // Display-only: histograms are NOT modified, no bins are filtered out.
        // Bins with h_cv == 0 are skipped from the max calculation because
        // the fractional shift is mathematically undefined there.
        double max_chi2_ndf = 0;
        double max_frac_shift = 0;
        int nbins_x = sh.h_cv->GetNbinsX();
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          double chi2 = sh.h_cv->Chi2Test(sh.h_var[i], "WW CHI2/NDF");
          if (std::isfinite(chi2) && chi2 > max_chi2_ndf) max_chi2_ndf = chi2;
          for (int b = 1; b <= nbins_x; ++b) {
            double cv_val = sh.h_cv->GetBinContent(b);
            if (cv_val == 0) continue;
            double var_val = sh.h_var[i]->GetBinContent(b);
            double frac = std::abs(var_val - cv_val) / std::abs(cv_val);
            if (frac > max_frac_shift) max_frac_shift = frac;
          }
        }

        // Draw stats box in upper-left of main pad
        TLatex stats;
        stats.SetNDC();
        stats.SetTextSize(0.055);
        stats.SetTextFont(42);
        stats.SetTextAlign(13); // top-left
        stats.DrawLatex(0.19, 0.92,
                         Form("#chi^{2}_{max}/ndf = %.2f", max_chi2_ndf));
        stats.DrawLatex(0.19, 0.85,
                         Form("max|#Delta|/CV = %.1f%%", max_frac_shift * 100));

        // Ratio pad
        c->cd();
        TPad *pratio = new TPad(Form("r%d", iv), "", x1, y_bot, x2, y_split);
        pratio->SetTopMargin(0.01); pratio->SetBottomMargin(0.30);
        pratio->SetLeftMargin(0.16); pratio->SetRightMargin(0.03);
        pratio->Draw(); pratio->cd();

        TH1D *rframe = (TH1D *)sh.h_cv->Clone(Form("rf%d", iv));
        rframe->Reset();
        rframe->SetMaximum(0.5); rframe->SetMinimum(-0.5);
        rframe->GetXaxis()->SetTitle(vd.label.c_str());
        rframe->GetXaxis()->SetTitleSize(0.14); rframe->GetXaxis()->SetTitleOffset(0.9);
        rframe->GetXaxis()->SetLabelSize(0.12);
        rframe->GetYaxis()->SetTitle("#frac{Var-CV}{CV}");
        rframe->GetYaxis()->SetTitleSize(0.12); rframe->GetYaxis()->SetTitleOffset(0.45);
        rframe->GetYaxis()->SetLabelSize(0.10);
        rframe->GetYaxis()->SetNdivisions(505);
        rframe->Draw("AXIS");

        TLine *zero = new TLine(rframe->GetXaxis()->GetXmin(), 0,
                                 rframe->GetXaxis()->GetXmax(), 0);
        zero->SetLineStyle(2); zero->SetLineColor(kGray + 1); zero->Draw();

        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          TH1D *rh = (TH1D *)sh.h_var[i]->Clone(Form("rh%d_%d", iv, i));
          rh->Add(sh.h_cv, -1.0);
          rh->Divide(sh.h_cv);
          int col_i, sty;
          GetVarColor(i, col_i, sty);
          rh->SetLineColor(col_i); rh->SetLineWidth(2); rh->SetLineStyle(sty);
          rh->Draw("HIST SAME");
        }
      }
      c->Print(pdfname.c_str());
      } // end ipage loop
    }
  }
  delete c;
}

// ===== ROOT output =========================================================
void WriteROOT(const HistMap &allhists) {
  std::string rootname = cliopts::output_base + ".root";
  TFile *fout = new TFile(rootname.c_str(), "RECREATE");

  for (auto &[ch, param_map] : allhists) {
    TDirectory *chdir = fout->mkdir(ch.c_str());
    for (auto &[par, var_map] : param_map) {
      TDirectory *pdir = chdir->mkdir(par.c_str());
      for (auto &[var, sh] : var_map) {
        TDirectory *vdir = pdir->mkdir(var.c_str());
        vdir->cd();
        sh.h_cv->Write("h_cv");
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          std::string hname = (i < (int)sh.tweakvals.size()) ?
            Form("h_var_%.3g", sh.tweakvals[i]) : Form("h_var_%d", i);
          sh.h_var[i]->Write(hname.c_str());
        }
      }
    }
  }

  // Save legend info as a TTree
  fout->cd();
  // (legend colors/styles are deterministic from GetVarColor, no need to store)

  fout->Close();
  std::cout << "ROOT file written to " << rootname << std::endl;
  delete fout;
}

// ===== main ================================================================
int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  if (cliopts::input_file.empty()) {
    std::cout << "[ERROR]: -i <input> is required." << std::endl;
    SayUsage(argv); return 1;
  }

  // Resolve variables
  std::vector<VarDef> vars;
  if (cliopts::variables.empty()) {
    for (auto &[key, vd] : kPredefinedVars) vars.push_back(vd);
  } else {
    for (auto &v : cliopts::variables) vars.push_back(ResolveVar(v));
  }

  // Detect input mode
  TFile *fin = TFile::Open(cliopts::input_file.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cout << "[ERROR]: Cannot open " << cliopts::input_file << std::endl;
    return 2;
  }

  bool flat_tree_mode = false;
  TTree *flat_tree = dynamic_cast<TTree *>(fin->Get(cliopts::treename.c_str()));
  TTree *meta_tree = dynamic_cast<TTree *>(fin->Get("tweak_metadata"));

  if (flat_tree && meta_tree) {
    TObjArray *branches = flat_tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i)
      if (std::string(branches->At(i)->GetName()).find("tweak_responses_") == 0)
        { flat_tree_mode = true; break; }
  }

  HistMap allhists;
  EventCounts counts;

  if (flat_tree_mode) {
    std::cout << "Detected flat tree mode." << std::endl;

    // Read param metadata
    std::vector<ParamMeta> params;
    {
      TObjString *name_obj = nullptr;
      int ntw = 0; double tweakvals[200];
      meta_tree->SetBranchAddress("name", &name_obj);
      meta_tree->SetBranchAddress("ntweaks", &ntw);
      meta_tree->SetBranchAddress("tweakvalues", tweakvals);
      for (Long64_t i = 0; i < meta_tree->GetEntries(); ++i) {
        meta_tree->GetEntry(i);
        ParamMeta pm;
        pm.fullname = name_obj->GetString().Data();
        pm.ntweaks = ntw;
        pm.tweakvalues.assign(tweakvals, tweakvals + ntw);
        if (ParamSelected(pm.fullname)) params.push_back(pm);
      }
      meta_tree->ResetBranchAddresses();
    }

    // Filter variables to those with existing branches (skip coslep etc.)
    std::vector<VarDef> valid_vars;
    for (auto &vd : vars) {
      bool is_computed = (vd.branch == "leading_proton_KE" || vd.branch == "coslep");
      if (is_computed) {
        if (vd.branch == "leading_proton_KE" && flat_tree->GetBranch("fsprotons_KE"))
          valid_vars.push_back(vd);
        else if (vd.branch == "coslep")
          std::cout << "[INFO]: coslep not available in flat tree mode, skipping." << std::endl;
        else
          std::cout << "[WARN]: " << vd.branch << " not available, skipping." << std::endl;
      } else if (flat_tree->GetBranch(vd.branch.c_str())) {
        valid_vars.push_back(vd);
      } else {
        std::cout << "[WARN]: Branch " << vd.branch << " not found, skipping." << std::endl;
      }
    }
    vars = valid_vars;

    FillFromFlatTree(flat_tree, meta_tree, vars, params, allhists, counts);
  } else {
    std::cout << "Detected GHEP mode." << std::endl;
    fin->Close(); delete fin; fin = nullptr;
    if (cliopts::fclname.empty()) {
      std::cout << "[ERROR]: GHEP mode requires -c <config.fcl>" << std::endl;
      return 3;
    }
    FillFromGHEP(cliopts::input_file, cliopts::fclname, vars, allhists, counts);
  }

  PruneEmpty(allhists);
  if (allhists.empty()) {
    std::cout << "[ERROR]: No histograms were filled." << std::endl;
    return 4;
  }

  ScaleToDiffXsec(allhists);

  // Build inclusive channels: total CC, total NC, total
  // by merging per-channel histograms
  std::set<std::string> param_names;
  for (auto &[ch, pm] : allhists)
    for (auto &[pn, _] : pm) param_names.insert(pn);

  auto AddChannel = [&](const std::string &newkey,
                         std::function<bool(const std::string &)> matcher) {
    for (auto &pn : param_names) {
      for (auto &vd : vars) {
        TH1D *sum_cv = nullptr;
        std::vector<TH1D *> sum_var;
        std::vector<double> tweakvals;

        for (auto &[ch, pm] : allhists) {
          if (!matcher(ch)) continue;
          auto pit = pm.find(pn);
          if (pit == pm.end()) continue;
          auto vit = pit->second.find(vd.branch);
          if (vit == pit->second.end()) continue;
          auto &sh = vit->second;

          if (!sum_cv) {
            std::string tag = newkey + "_" + pn + "_" + vd.branch;
            sum_cv = (TH1D *)sh.h_cv->Clone(("h_cv_" + tag).c_str());
            tweakvals = sh.tweakvals;
            for (int i = 0; i < (int)sh.h_var.size(); ++i)
              sum_var.push_back((TH1D *)sh.h_var[i]->Clone(Form("h_var%d_%s", i, tag.c_str())));
          } else {
            sum_cv->Add(sh.h_cv);
            for (int i = 0; i < (int)sh.h_var.size() && i < (int)sum_var.size(); ++i)
              sum_var[i]->Add(sh.h_var[i]);
          }
        }
        if (sum_cv && sum_cv->GetEntries() > 0) {
          SystHists sh;
          sh.h_cv = sum_cv;
          sh.h_var = sum_var;
          sh.tweakvals = tweakvals;
          allhists[newkey][pn][vd.branch] = sh;
        }
      }
    }
  };

  AddChannel("Total", [](const std::string &) { return true; });
  AddChannel("Total_CC", [](const std::string &ch) { return ch.find("CC_") == 0; });
  AddChannel("Total_NC", [](const std::string &ch) {
    return ch.find("NC_") == 0 && ch.find("Total") == std::string::npos;
  });

  if (cliopts::make_root) WriteROOT(allhists);

  if (cliopts::make_pdf) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    std::string pdfname = cliopts::output_base + ".pdf";
    TCanvas *c = new TCanvas("c", "", 1400, 1000);
    c->Print((pdfname + "[").c_str());

    // Summary page
    MakeSummaryPage(c, counts);
    c->Print(pdfname.c_str());

    // Now delegate to MakePlots for the actual histogram pages
    // (MakePlots writes pages to an already-open PDF)
    MakePlots(allhists, vars);

    c->Print((pdfname + "]").c_str());
    delete c;
  }

  if (fin) { fin->Close(); delete fin; }
  return 0;
}
