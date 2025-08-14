/*******************************************************************************
 * ValenciaMECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "ValenciaMECq0q3InterpWeighting_tool.hh"
#include <fhiclcpp/ParameterSet.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TString.h>  
#include <iostream>
#include <cmath>

// SystematicsTools includes
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

// GENIE includes for particle access
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace nusyst;

// ---------------------------------------------------------------------------
ValenciaMECq0q3InterpWeighting::ValenciaMECq0q3InterpWeighting(
    const fhicl::ParameterSet& p)
  : IGENIESystProvider_tool(p)  // Call base class constructor
{
  // ---------------------------------------------------------------------
  // Energy grid: support either new key "EnergyGrid" or legacy "Energies"
  // ---------------------------------------------------------------------
  if (p.has_key("EnergyGrid")) {
    fEgrid = p.get<std::vector<double>>("EnergyGrid");
  } else if (p.has_key("Energies")) {
    fEgrid = p.get<std::vector<double>>("Energies");
  } else {
    throw std::runtime_error("ValenciaMECq0q3InterpWeighting: Need EnergyGrid (or legacy Energies) in config");
  }
  if (fEgrid.empty()) throw std::runtime_error("ValenciaMECq0q3InterpWeighting: Provided energy grid is empty");

  // ---------------------------------------------------------------------
  // Weight limits: support table WeightLimits: { min: X, max: Y } OR list
  // WeightLimits : [ X , Y ]
  // ---------------------------------------------------------------------
  fWmin = 0.0; fWmax = 5.0; // defaults
  if (p.has_key("WeightLimits.min") || p.has_key("WeightLimits.max")) {
    fWmin = p.get<double>("WeightLimits.min", fWmin);
    fWmax = p.get<double>("WeightLimits.max", fWmax);
  } else if (p.has_key("WeightLimits")) {
    try {
      // Try as table first (above would have caught min/max), so interpret as sequence
      std::vector<double> wl = p.get<std::vector<double>>("WeightLimits");
      if (wl.size() >= 2) { fWmin = wl.front(); fWmax = wl.back(); }
    } catch (...) {
      // fall back silently, keep defaults
    }
  }
  if (fWmin <= 0.0) fWmin = 0.000; // guard against zero for reciprocal ops later if any
  if (fWmax < fWmin) std::swap(fWmin, fWmax);

  // ---------------------------------------------------------------------
  // Histogram sourcing strategies:
  //  (A) Single WeightFile containing all h_weights_map_[np/nn]_X.YGeV
  //  (B) Arrays np_files / nn_files each parallel to energy grid
  // ---------------------------------------------------------------------
  bool haveSingle = p.has_key("WeightFile");
  bool haveArrays = p.has_key("np_files") && p.has_key("nn_files");
  if (!haveSingle && !haveArrays) {
    throw std::runtime_error("ValenciaMECq0q3InterpWeighting: Need either WeightFile or (np_files & nn_files)");
  }
  if (haveSingle && haveArrays) {
    // Ambiguous, prefer explicit WeightFile
    haveArrays = false;
  }

  if (haveSingle) {
    const std::string fname = p.get<std::string>("WeightFile");
  fToolOptions.put("WeightFile", fname);
    TFile f(fname.c_str(), "READ");
    if (!f.IsOpen()) throw std::runtime_error("Cannot open WeightFile " + fname);
    for (Topo topo : {Topo::np, Topo::nn}) {
      const char* tag = (topo == Topo::np ? "np" : "nn");
      auto& v = fCalcs[topo]; v.reserve(fEgrid.size());
      for (double E : fEgrid) {
        const std::string hname = Form("h_weights_map_%s_%0.1fGeV", tag, E);
        if (TH2D* h = dynamic_cast<TH2D*>(f.Get(hname.c_str()))) {
          v.emplace_back(std::make_unique<ValenciaMECq0q3ResponseCalc>(h, fWmin, fWmax));
        } else {
          throw std::runtime_error("Missing histogram " + hname + " in " + fname);
        }
      }
    }
  } else { // haveArrays
    auto np_files = p.get<std::vector<std::string>>("np_files");
    auto nn_files = p.get<std::vector<std::string>>("nn_files");
  fToolOptions.put("np_files", np_files);
  fToolOptions.put("nn_files", nn_files);
    if (np_files.size() != fEgrid.size() || nn_files.size() != fEgrid.size()) {
      throw std::runtime_error("ValenciaMECq0q3InterpWeighting: np_files/nn_files size must match energy grid size");
    }
    // Load each file individually
    auto load = [&](const std::vector<std::string>& files, Topo topo) {
      const char* tag = (topo == Topo::np ? "np" : "nn");
      auto& v = fCalcs[topo]; v.reserve(fEgrid.size());
      for (size_t i = 0; i < files.size(); ++i) {
        const double E = fEgrid[i];
        const std::string& fname = files[i];
        TFile f(fname.c_str(), "READ");
        if (!f.IsOpen()) throw std::runtime_error("Cannot open file " + fname);
        const std::string hname = Form("h_weights_map_%s_%0.1fGeV", tag, E);
        TH2D* h = dynamic_cast<TH2D*>(f.Get(hname.c_str()));
        if (!h) {
          // Try a fallback: maybe histogram does not include energy in name, just one per file
          // Accept first TH2D in file if expected name missing
          TIter nextkey(f.GetListOfKeys());
          TH2D* firstH = nullptr;
          while (TKey* key = (TKey*)nextkey()) {
            TObject* obj = key->ReadObj();
            if ((firstH = dynamic_cast<TH2D*>(obj))) break;
          }
          if (!firstH) {
            throw std::runtime_error("Missing histogram " + hname + " (and no fallback TH2D) in " + fname);
          }
          h = firstH; // use fallback
        }
        v.emplace_back(std::make_unique<ValenciaMECq0q3ResponseCalc>(h, fWmin, fWmax));
      }
    };
    load(np_files, Topo::np);
    load(nn_files, Topo::nn);
  }

  // Store common options
  fToolOptions.put("Energies", fEgrid);
  // Preserve original WeightLimits form for output (vector)
  fToolOptions.put("WeightLimits", std::vector<double>{fWmin, fWmax});
}

// ---------------------------------------------------------------------------
systtools::SystMetaData ValenciaMECq0q3InterpWeighting::BuildSystMetaData(
    fhicl::ParameterSet const &ps, systtools::paramId_t firstId) {
  
  systtools::SystMetaData smd;
  
  // Define a single parameter for Valencia MEC reweighting
  systtools::SystParamHeader phdr;
  // Prefer new key ValenciaMECResponse; fall back to legacy ValenciaMEC_q0q3Interp
  if (systtools::ParseFhiclToolConfigurationParameter(ps, "ValenciaMECResponse", phdr, firstId) ||
      systtools::ParseFhiclToolConfigurationParameter(ps, "ValenciaMEC_q0q3Interp", phdr, firstId)) {
    phdr.systParamId = firstId++;
    // If no prettyName provided, set a reasonable default
    if (phdr.prettyName.empty()) phdr.prettyName = "ValenciaMEC q0q3 Interp";
    smd.push_back(phdr);
  } else {
    // Optional parameter: if absent, tool will yield unity weights (no variations)
    std::cerr << "[ValenciaMECq0q3InterpWeighting] Warning: No ValenciaMECResponse (or legacy"
              << " ValenciaMEC_q0q3Interp) parameter header provided; returning unity weights." << std::endl;
  }
  
  return smd;
}

// ---------------------------------------------------------------------------
bool ValenciaMECq0q3InterpWeighting::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  
  // Store the metadata for later use
  systtools::SystMetaData const &md = GetSystMetaData();
  
  // This class handles its own response calculation in the constructor
  // and GetEventResponse method, so just return true
  return true;
}

// ---------------------------------------------------------------------------
systtools::event_unit_response_t ValenciaMECq0q3InterpWeighting::GetEventResponse(
    genie::EventRecord const& ev)
{
  // ---  classify topology -------------------------------------------------
  const Topo topo = ClassifyEvent(ev);
  if (topo == Topo::unknown) {
    // Return response for this parameter ID with default response
    systtools::event_unit_response_t response;
    if (!this->GetSystMetaData().empty()) {
      systtools::ParamResponses pr;
      pr.pid = this->GetSystMetaData()[0].systParamId;
      pr.responses = {1.0, 1.0, 1.0};  // down, central, up
      response.push_back(pr);
    }
    return response;
  }

  // ---  compute kinematics -------------------------------------------------
  double q0, q3, Enu;
  ComputeQ0Q3(ev, q0, q3, Enu);

  // ---  find bracketing energy indices ------------------------------------
  auto it_high = std::lower_bound(fEgrid.begin(), fEgrid.end(), Enu);
  size_t i_high = (it_high == fEgrid.end() ? fEgrid.size() - 1
                                           : std::distance(fEgrid.begin(), it_high));
  size_t i_low  = (i_high == 0 ? 0 : i_high - 1);
  const double Elo = fEgrid[i_low];
  const double Ehi = fEgrid[i_high];
  const double t   = (Ehi > Elo) ? (Enu - Elo) / (Ehi - Elo) : 0.0;

  // ---  bilinear map lookup + linear blend --------------------------------
  const auto& vec = fCalcs.at(topo);
  double w_lo = vec[i_low ]->GetCentralWeight(q0, q3);
  double w_hi = vec[i_high]->GetCentralWeight(q0, q3);
  double w_cv = (1.0 - t) * w_lo + t * w_hi; // interpolated central value
  w_cv = std::clamp(w_cv, fWmin, fWmax);

  systtools::event_unit_response_t response;
  if (!this->GetSystMetaData().empty()) {
    const auto &hdr = this->GetSystMetaData()[0];
    systtools::ParamResponses pr;
    pr.pid = hdr.systParamId;
    pr.responses.reserve(hdr.paramVariations.size());
    double delta = w_cv - 1.0; // deviation from unity
    if (std::abs(delta) < 1e-12) delta = 0.0; // guard tiny noise
    for (double d : hdr.paramVariations) {
      double w = w_cv + d * delta; // linear scaling around unity keeping d=0 => w_cv
      w = std::clamp(w, fWmin, fWmax);
      pr.responses.push_back(w);
    }
    response.push_back(pr);
  }
  return response;
}

/*******************************************************************************
 *  Additional helper implementations for ValenciaMECq0q3InterpWeighting_tool.cc
 *  This file contains the implementation of methods for extracting kinematics
 *  and classifying event topologies.
 ******************************************************************************/

// ---------------------------------------------------------------------------
//  Extract (q0, q3, Enu) from the GENIE event.
// ---------------------------------------------------------------------------
void ValenciaMECq0q3InterpWeighting::ComputeQ0Q3(genie::EventRecord const& ev,
                                                 double& q0, double& q3,
                                                 double& Enu)
{
    const TLorentzVector p4nu = *ev.Probe()->P4();
    Enu = p4nu.E();

    const TLorentzVector p4lep = *ev.FinalStatePrimaryLepton()->P4();
    TLorentzVector qv = p4nu - p4lep;   // four‑vector transfer
    q0 = qv.E();
    q3 = qv.P();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//  Looks at nucleon-target particles (status kIStNucleonTarget).
// ---------------------------------------------------------------------------
ValenciaMECq0q3InterpWeighting::Topo
ValenciaMECq0q3InterpWeighting::ClassifyEvent(genie::EventRecord const& ev)
{
  for(int i=0;i<ev.GetEntries();++i){
    const auto* p=ev.Particle(i);
    if(!p) continue;
    if(p->Status()!=genie::kIStNucleonTarget) continue;
    if(p->Pdg()==genie::kPdgClusterNN) return Topo::nn; // 2n cluster
    if(p->Pdg()==genie::kPdgClusterNP) return Topo::np; // 1n+1p cluster
  }
  return Topo::unknown;
}
