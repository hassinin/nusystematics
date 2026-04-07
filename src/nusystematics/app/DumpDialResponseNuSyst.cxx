// DumpDialResponseNuSyst.cxx
// For each configured dial, sample N events per interaction channel and
// plot the weight-vs-parameter-value curve per event. Reveals the shape of
// the dial response (linearity, monotonicity, asymmetry, event dependence).
//
// Uses IGENIESystProvider_tool::GetEventResponse(evt, paramVals) to force
// arbitrary grid evaluation via OverrideVariations.

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystParamHeader.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/response_helper.hh"

#include "fhiclcpp/ParameterSet.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace systtools;
using namespace nusyst;

// ===== CLI =================================================================
namespace cliopts {
std::string fclname;
std::string input_file;
std::string output_base = "dial_response";
std::string genie_branch_name = "gmcrec";
std::vector<std::string> parameters;
size_t NMax = std::numeric_limits<size_t>::max();
int n_per_channel = 5;
int nsteps = 0; // 0 = use paramVariations as-is
bool make_pdf = true;
bool make_root = true;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Required:\n"
    "    -c <config.fcl>      FHiCL config file\n"
    "    -i <ghep.root>       GENIE ghep input\n\n"
    "  Optional:\n"
    "    -p <par1,par2,...>    Filter dials by name substring\n"
    "    -n <N>                Events per channel (default 5)\n"
    "    --nsteps <N>          Resample paramVariations range with N points\n"
    "                          (default: use paramVariations as-is)\n"
    "    -o <basename>         Output basename (default: dial_response)\n"
    "    -b <branch>           NtpMCEventRecord branch (default: gmcrec)\n"
    "    -N <NMax>             Max events to scan for sampling\n"
    "    --no-pdf / --no-root  Skip output\n"
    "    -?|--help             Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "-c") cliopts::fclname = argv[++opt];
    else if (s == "-i") cliopts::input_file = argv[++opt];
    else if (s == "-o") cliopts::output_base = argv[++opt];
    else if (s == "-b") cliopts::genie_branch_name = argv[++opt];
    else if (s == "-n") cliopts::n_per_channel = str2T<int>(argv[++opt]);
    else if (s == "-N") cliopts::NMax = str2T<size_t>(argv[++opt]);
    else if (s == "--nsteps") cliopts::nsteps = str2T<int>(argv[++opt]);
    else if (s == "--no-pdf") cliopts::make_pdf = false;
    else if (s == "--no-root") cliopts::make_root = false;
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

// ===== Channel classification (copied from PlotSystVariationsNuSyst) =======
std::string NuName(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "numu" : "numubar";
    case 12: return pdg > 0 ? "nue"  : "nuebar";
    case 16: return pdg > 0 ? "nutau" : "nutaubar";
    default: return Form("nu%d", pdg);
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

std::string InteractionType(const genie::EventRecord &ev) {
  auto const &proc = ev.Summary()->ProcInfo();
  if (proc.IsQuasiElastic())    return "QE";
  if (proc.IsMEC())             return "MEC";
  if (proc.IsResonant())        return "RES";
  if (proc.IsDeepInelastic())   return "DIS";
  return "Other";
}

std::string MakeChannelKey(const genie::EventRecord &ev) {
  std::string cc = ev.Summary()->ProcInfo().IsWeakCC() ? "CC" : "NC";
  int nu_pdg = ev.Probe() ? ev.Probe()->Pdg() : 0;
  int tgt_A = ev.TargetNucleus() ? ev.TargetNucleus()->A() : 1;
  return cc + "_" + NuName(nu_pdg) + "_" + TargetName(tgt_A) + "_" + InteractionType(ev);
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

// ===== Data containers =====================================================
struct DialInfo {
  paramId_t pid;
  size_t provider_idx;
  std::string provider_name;
  std::string name;
  std::string fullname;
  double central;
  bool is_correction;
  bool is_natural;
  std::vector<double> grid; // x-values (sigma or natural) where we evaluate
};

struct EventRef {
  Long64_t index;
  std::string channel;
};

// Per (dial, event) → vector of weights aligned with dial.grid
using ResponseMap = std::map<std::string, std::map<std::string, std::vector<std::pair<Long64_t, std::vector<double>>>>>;
// ResponseMap[dial_fullname][channel] -> list of (event_index, weights)

// ===== Main ================================================================
int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  if (cliopts::fclname.empty() || cliopts::input_file.empty()) {
    std::cout << "[ERROR]: -c and -i are required." << std::endl;
    SayUsage(argv);
    return 1;
  }

  // Load providers
  response_helper phh(cliopts::fclname);

  // Collect dial info
  std::vector<DialInfo> dials;
  for (paramId_t pid : phh.GetParameters()) {
    SystParamHeader const &hdr = phh.GetHeader(pid);
    if (hdr.isResponselessParam) continue;
    if (!ParamSelected(hdr.prettyName)) continue;

    // Find provider
    size_t prov_idx = 0;
    std::string provname = "?";
    for (size_t si = 0; si < phh.GetSystProvider().size(); ++si) {
      if (phh.GetSystProvider()[si]->ParamIsHandled(pid)) {
        prov_idx = si;
        provname = phh.GetSystProvider()[si]->GetFullyQualifiedName();
        break;
      }
    }

    DialInfo di;
    di.pid = pid;
    di.provider_idx = prov_idx;
    di.provider_name = provname;
    di.name = hdr.prettyName;
    di.fullname = provname + "_" + hdr.prettyName;
    di.central = hdr.centralParamValue;
    di.is_correction = hdr.isCorrection;
    di.is_natural = hdr.unitsAreNatural;

    if (hdr.isCorrection || hdr.paramVariations.empty()) {
      std::cout << "[INFO]: Skipping " << di.fullname
                << " (correction or no variations)." << std::endl;
      continue;
    }

    // Build evaluation grid
    if (cliopts::nsteps > 0) {
      double vmin = *std::min_element(hdr.paramVariations.begin(), hdr.paramVariations.end());
      double vmax = *std::max_element(hdr.paramVariations.begin(), hdr.paramVariations.end());
      di.grid.reserve(cliopts::nsteps);
      for (int i = 0; i < cliopts::nsteps; ++i) {
        double t = (cliopts::nsteps == 1) ? 0.5 : double(i) / (cliopts::nsteps - 1);
        di.grid.push_back(vmin + t * (vmax - vmin));
      }
    } else {
      di.grid = hdr.paramVariations;
    }

    dials.push_back(std::move(di));
  }

  if (dials.empty()) {
    std::cout << "[ERROR]: No dials to process." << std::endl;
    return 2;
  }

  std::cout << "Configured " << dials.size() << " dials." << std::endl;

  // Open ghep input
  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(cliopts::input_file.c_str())) {
    std::cout << "[ERROR]: Failed to add " << cliopts::input_file << std::endl;
    return 3;
  }
  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  gevs->SetBranchAddress(cliopts::genie_branch_name.c_str(), &GenieNtpl);

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");

  // Pass 1: scan events to pick N per channel
  Long64_t nentries = std::min((Long64_t)cliopts::NMax, gevs->GetEntries());
  std::map<std::string, std::vector<Long64_t>> channel_events;

  std::cout << "Scanning " << nentries << " events to select samples..." << std::endl;
  for (Long64_t ev = 0; ev < nentries; ++ev) {
    gevs->GetEntry(ev);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;
    std::string chkey = MakeChannelKey(GenieGHep);
    auto &vec = channel_events[chkey];
    if ((int)vec.size() < cliopts::n_per_channel) {
      vec.push_back(ev);
    }
    GenieNtpl->Clear();
  }

  size_t total_selected = 0;
  for (auto &[ch, evs] : channel_events) total_selected += evs.size();
  std::cout << "Selected " << total_selected << " events across "
            << channel_events.size() << " channels." << std::endl;

  // Pass 2: for each selected event, compute dial response curves
  // Storage: responses[dial_fullname][channel] -> list of (event_index, weights)
  ResponseMap responses;

  std::cout << "Computing dial responses..." << std::endl;
  size_t processed = 0;
  for (auto &[chkey, event_list] : channel_events) {
    for (Long64_t ev_idx : event_list) {
      gevs->GetEntry(ev_idx);
      genie::EventRecord const &GenieGHep = *GenieNtpl->event;

      for (auto &di : dials) {
        IGENIESystProvider_tool *prov = phh.GetSystProvider()[di.provider_idx].get();

        // Build paramVals for this single dial with its grid
        std::vector<std::pair<paramId_t, std::vector<double>>> paramVals = {
            {di.pid, di.grid}};

        std::vector<double> weights(di.grid.size(), 1.0);
        try {
          auto resp = prov->GetEventResponse(GenieGHep, paramVals);
          // Find our pid in the response
          for (auto const &pr : resp) {
            if (pr.pid == di.pid) {
              for (size_t i = 0; i < pr.responses.size() && i < weights.size(); ++i)
                weights[i] = pr.responses[i];
              break;
            }
          }
        } catch (std::exception &e) {
          // leave weights at 1
        }

        responses[di.fullname][chkey].push_back({ev_idx, std::move(weights)});
      }

      processed++;
      GenieNtpl->Clear();
    }
  }
  std::cout << "Processed " << processed << " events." << std::endl;

  // Write ROOT output
  if (cliopts::make_root) {
    std::string rootname = cliopts::output_base + ".root";
    TFile *fout = new TFile(rootname.c_str(), "RECREATE");
    for (auto &di : dials) {
      TDirectory *ddir = fout->mkdir(di.fullname.c_str());
      auto dit = responses.find(di.fullname);
      if (dit == responses.end()) continue;
      for (auto &[chkey, ev_list] : dit->second) {
        TDirectory *cdir = ddir->mkdir(chkey.c_str());
        cdir->cd();
        for (size_t i = 0; i < ev_list.size(); ++i) {
          auto &[ev_idx, w] = ev_list[i];
          TGraph *g = new TGraph((int)di.grid.size());
          for (size_t k = 0; k < di.grid.size(); ++k)
            g->SetPoint(k, di.grid[k], w[k]);
          g->SetName(Form("event_%lld", ev_idx));
          g->SetTitle(Form("%s | %s | event %lld;%s;weight",
                            di.fullname.c_str(), chkey.c_str(), ev_idx,
                            di.is_natural ? "param value (natural)" : "param value (sigma)"));
          g->Write();
        }
      }
    }
    fout->Close();
    delete fout;
    std::cout << "ROOT file written to " << rootname << std::endl;
  }

  // Write PDF
  if (cliopts::make_pdf) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::string pdfname = cliopts::output_base + ".pdf";
    TCanvas *c = new TCanvas("c", "", 1200, 900);
    c->Print((pdfname + "[").c_str());

    // Summary page: list sampled events per channel
    c->Clear();
    TPad *sumpad = new TPad("sumpad", "", 0, 0, 1, 1);
    sumpad->SetFillColor(kWhite);
    sumpad->Draw(); sumpad->cd();
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.05);
    tex.SetTextFont(62);
    tex.DrawLatex(0.05, 0.93, "DumpDialResponseNuSyst Summary");
    tex.SetTextFont(42);
    tex.SetTextSize(0.028);
    double y = 0.86;
    tex.DrawLatex(0.05, y, Form("Config: %s", cliopts::fclname.c_str())); y -= 0.03;
    tex.DrawLatex(0.05, y, Form("Input:  %s", cliopts::input_file.c_str())); y -= 0.03;
    tex.DrawLatex(0.05, y, Form("Dials:  %zu   Channels: %zu   Events/channel: %d",
                                 dials.size(), channel_events.size(), cliopts::n_per_channel));
    y -= 0.04;
    tex.SetTextFont(62);
    tex.DrawLatex(0.05, y, "Sampled events per channel:"); y -= 0.03;
    tex.SetTextFont(42);
    for (auto &[chkey, evs] : channel_events) {
      std::ostringstream ss;
      ss << chkey << ":  ";
      for (size_t i = 0; i < evs.size(); ++i) {
        if (i > 0) ss << ", ";
        ss << evs[i];
      }
      tex.DrawLatex(0.07, y, ss.str().c_str());
      y -= 0.025;
      if (y < 0.05) break;
    }
    c->Print(pdfname.c_str());

    // One page per (dial, channel)
    for (auto &di : dials) {
      auto dit = responses.find(di.fullname);
      if (dit == responses.end()) continue;

      for (auto &[chkey, ev_list] : dit->second) {
        if (ev_list.empty()) continue;

        c->Clear();
        c->cd();

        // Title pad
        TPad *tp = new TPad("tp", "", 0, 0.93, 1, 1.0);
        tp->SetFillColor(kWhite);
        tp->Draw(); tp->cd();
        TLatex title;
        title.SetNDC(); title.SetTextSize(0.45); title.SetTextAlign(22); title.SetTextFont(62);
        std::string ptitle = di.fullname + "  |  " + chkey;
        title.DrawLatex(0.5, 0.5, ptitle.c_str());

        // Main pad
        c->cd();
        TPad *mp = new TPad("mp", "", 0.08, 0.04, 0.98, 0.92);
        mp->SetFillColor(kWhite);
        mp->SetLeftMargin(0.13);
        mp->SetBottomMargin(0.13);
        mp->SetTopMargin(0.05);
        mp->SetRightMargin(0.05);
        mp->Draw(); mp->cd();

        // Find weight range
        double wmin = 1, wmax = 1;
        for (auto &[idx, w] : ev_list) {
          for (double v : w) {
            if (std::isfinite(v)) {
              wmin = std::min(wmin, v);
              wmax = std::max(wmax, v);
            }
          }
        }
        double pad = 0.1 * (wmax - wmin);
        if (pad == 0) pad = 0.1;
        double gmin = *std::min_element(di.grid.begin(), di.grid.end());
        double gmax = *std::max_element(di.grid.begin(), di.grid.end());

        // Frame histogram
        TH1D *frame = new TH1D(Form("frame_%s_%s", di.fullname.c_str(), chkey.c_str()),
                                 "", 100, gmin, gmax);
        frame->SetMinimum(wmin - pad);
        frame->SetMaximum(wmax + pad);
        std::string xaxis = di.is_natural ? "parameter value (natural units)"
                                           : "parameter value (#sigma)";
        frame->GetXaxis()->SetTitle(xaxis.c_str());
        frame->GetYaxis()->SetTitle("weight");
        frame->GetXaxis()->SetTitleSize(0.045);
        frame->GetYaxis()->SetTitleSize(0.045);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->Draw("AXIS");

        // Unity line (weight=1)
        TLine *unity = new TLine(gmin, 1, gmax, 1);
        unity->SetLineStyle(2); unity->SetLineColor(kGray + 1);
        unity->Draw();

        // CV vertical line
        if (di.central != kDefaultDouble && di.central >= gmin && di.central <= gmax) {
          TLine *cvline = new TLine(di.central, wmin - pad, di.central, wmax + pad);
          cvline->SetLineStyle(3); cvline->SetLineColor(kGray + 2);
          cvline->Draw();
        }

        // One TGraph per event
        TLegend *leg = new TLegend(0.65, 0.65, 0.94, 0.93);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextSize(0.032);
        leg->SetHeader("event index");

        for (size_t i = 0; i < ev_list.size(); ++i) {
          auto &[ev_idx, w] = ev_list[i];
          TGraph *g = new TGraph((int)di.grid.size());
          for (size_t k = 0; k < di.grid.size(); ++k)
            g->SetPoint(k, di.grid[k], w[k]);
          int col, sty;
          GetVarColor((int)i, col, sty);
          g->SetLineColor(col);
          g->SetLineWidth(2);
          g->SetLineStyle(sty);
          g->SetMarkerColor(col);
          g->SetMarkerStyle(20 + (int)(i % 8));
          g->SetMarkerSize(0.9);
          g->Draw("LP SAME");
          leg->AddEntry(g, Form("%lld", ev_idx), "lp");
        }
        leg->Draw();

        c->Print(pdfname.c_str());
      }
    }
    c->Print((pdfname + "]").c_str());
    std::cout << "PDF written to " << pdfname << std::endl;
    delete c;
  }

  delete gevs;
  return 0;
}
