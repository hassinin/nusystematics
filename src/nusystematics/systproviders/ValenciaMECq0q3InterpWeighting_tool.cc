/*******************************************************************************
 * ValenciaMECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "ValenciaMECq0q3InterpWeighting_tool.hh"
#include <fhiclcpp/ParameterSet.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TString.h>  // For Form function
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
  : IGENIESystProvider_tool(p),  // Call base class constructor
    fEgrid(p.get<std::vector<double>>("EnergyGrid")),
    fWmin(p.get<double>("WeightLimits.min", 0.1)),
    fWmax(p.get<double>("WeightLimits.max", 5.0))
{
  const std::string fname = p.get<std::string>("WeightFile");
  TFile f(fname.c_str(), "READ");
  if (!f.IsOpen()) throw std::runtime_error("Cannot open " + fname);

  for (Topo topo : {Topo::np, Topo::nn}) {
    const char* tag = (topo == Topo::np ? "np" : "nn");
    auto& v = fCalcs[topo]; v.reserve(fEgrid.size());
    for (double E : fEgrid) {
      const std::string hname = Form("h_weights_map_%s_%0.1fGeV", tag, E);
      if (TH2D* h = dynamic_cast<TH2D*>(f.Get(hname.c_str()))) {
        v.emplace_back(std::make_unique<ValenciaMECq0q3ResponseCalc>(h, fWmin, fWmax));
      } else {
        throw std::runtime_error("Missing histogram " + hname);
      }
    }
  }
}

// ---------------------------------------------------------------------------
systtools::SystMetaData ValenciaMECq0q3InterpWeighting::BuildSystMetaData(
    fhicl::ParameterSet const &ps, systtools::paramId_t firstId) {
  
  systtools::SystMetaData smd;
  
  // Define a single parameter for Valencia MEC reweighting
  systtools::SystParamHeader phdr;
  if (systtools::ParseFhiclToolConfigurationParameter(ps, "ValenciaMECResponse", phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
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
  double w_cv = (1.0 - t) * w_lo + t * w_hi;
  double w_up = std::clamp(2.0 - w_cv, fWmin, fWmax);
  double w_dn = w_up; // symmetric

  // Package into the expected response format
  systtools::event_unit_response_t response;
  if (!this->GetSystMetaData().empty()) {
    systtools::ParamResponses pr;
    pr.pid = this->GetSystMetaData()[0].systParamId;
    pr.responses = {w_dn, w_cv, w_up};
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
