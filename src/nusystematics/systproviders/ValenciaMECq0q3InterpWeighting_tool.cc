//======================================================================
//  ValenciaMECq0q3InterpWeighting_tool.cc   (re-write)
//======================================================================
#include "ValenciaMECq0q3InterpWeighting_tool.hh"
#include "nusystematics/responsecalculators/ValenciaMECq0q3ResponseCalc.hh"

#include "systematicstools/interface/EventResponse_product.hh"
#include "nusystematics/utility/exceptions.hh"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "fhiclcpp/ParameterSet.h"
#include "Framework/GHEP/GHepParticle.h"
#include "TLorentzVector.h"
#include "Framework/ParticleData/PDGCodes.h"
#include <algorithm>
#include <numeric>

using namespace systtools;   // only inside this translation unit

namespace nusyst
{
//--------------------------------------------------------------------
// ctor
ValenciaMECq0q3InterpWeighting::
ValenciaMECq0q3InterpWeighting(fhicl::ParameterSet const &ps)
  : IGENIESystProvider_tool(ps) {}

//--------------------------------------------------------------------
// metadata
SystMetaData
ValenciaMECq0q3InterpWeighting::BuildSystMetaData(fhicl::ParameterSet const &cfg,
                                                  paramId_t firstId)
{
  SystMetaData md;
  SystParamHeader h;
  h.systParamId       = firstId++;
  h.prettyName        = "ValenciaMEC_q0q3Interp";
  h.centralParamValue = 0.0;
  h.paramVariations   = {-1.,0.,1.};
  
  md.push_back(h);

  // copy only simple values / vectors – no nested tables
  if (cfg.has_key("Energies"))
    tool_options.put("Energies", cfg.get<std::vector<double>>("Energies"));

  if (cfg.has_key("WeightLimits"))
    tool_options.put("WeightLimits", cfg.get<std::vector<double>>("WeightLimits"));

  if (cfg.has_key("np_files"))
    tool_options.put("np_files", cfg.get<std::vector<std::string>>("np_files"));

  if (cfg.has_key("nn_files"))
    tool_options.put("nn_files", cfg.get<std::vector<std::string>>("nn_files"));
  return md;
}

//--------------------------------------------------------------------
// set-up (load histograms & create calculators)
bool
ValenciaMECq0q3InterpWeighting::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_opts)
{
  // Read directly from tool_opts instead of looking for nested parameter
  fhicl::ParameterSet const &pset = tool_opts;

  // energies ---------------------------------------------------------
  fEnergies = pset.get<std::vector<double>>("Energies");
  if (fEnergies.empty())
    throw invalid_ToolConfigurationFHiCL() << "[ValenciaMEC] missing Energies";
  std::sort(fEnergies.begin(), fEnergies.end());

  // helper to load one topology -------------------------------------
  auto load = [&](Topology topo, std::string const &key)
  {
    if (!pset.has_key(key)) return;
    auto const files = pset.get<std::vector<std::string>>(key);
    fCalcs[topo].resize(fEnergies.size());

    for (std::size_t i = 0; i < fEnergies.size(); ++i) {
      if (i >= files.size()) continue;
      TFile f(files[i].c_str(), "READ");
      if (!f.IsOpen()) continue;
      TH2 *raw = nullptr;
      f.GetObject("h_weights_map", raw);
      if (!raw) continue;
      std::unique_ptr<TH2> h(static_cast<TH2*>(raw->Clone()));
      h->SetDirectory(nullptr);
      fCalcs[topo][i] =
        std::make_unique<ValenciaMECq0q3ResponseCalc>(h.release());
    }
  };
  load(Topology::np, "np_files");
  load(Topology::nn, "nn_files");

  // optional limits --------------------------------------------------
  if (pset.has_key("WeightLimits")) {
    auto v = pset.get<std::vector<double>>("WeightLimits");
    if (v.size() == 2) fWeightLimits = {v[0], v[1]};
  }
  return true;
}

//--------------------------------------------------------------------
// event helpers
Topology
ValenciaMECq0q3InterpWeighting::classifyEvent(genie::EventRecord const &ev) const
{
  for (int i = 0; i < ev.GetEntries(); ++i) {
  auto const *p = ev.Particle(i); // GHepParticle*
  if (!p) continue;
  if (p->Status() == genie::kIStNucleonTarget) {
    if (p->Pdg() == genie::kPdgClusterNN) return Topology::nn;
    if (p->Pdg() == genie::kPdgClusterNP) return Topology::np;
  }
}
return Topology::unknown;
}
//--------------------------------------------------------------------
void
ValenciaMECq0q3InterpWeighting::computeQ0Q3(genie::EventRecord const &ev,
                                            double &q0,double &q3) const
{
  auto *in  = ev.Probe();
  auto *out = ev.FinalStatePrimaryLepton();
  if (!in || !out) { q0 = q3 = 0.; return; }
  TLorentzVector t = *(in->P4()) - *(out->P4());
  q0 = t.E();
  q3 = t.Vect().Mag();
}
//--------------------------------------------------------------------
double
ValenciaMECq0q3InterpWeighting::getInterpolatedWeight(double E,double q0,double q3,
                                                      Topology topo) const
{
  auto itTopo = fCalcs.find(topo);
  if (itTopo == fCalcs.end()) return 1.0;

  auto const &set = itTopo->second;
  if (set.empty()) return 1.0;

  // below / above grid
  if (E <= fEnergies.front())
    return set.front()
             ? set.front()->GetVariation(0,set.front()->GetBin({q0,q3})) : 1.0;
  if (E >= fEnergies.back())
    return set.back()
             ? set.back()->GetVariation(0,set.back()->GetBin({q0,q3})) : 1.0;

  // interpolate
  auto itHi = std::upper_bound(fEnergies.begin(), fEnergies.end(), E);
  std::size_t iHi = itHi - fEnergies.begin();
  std::size_t iLo = iHi - 1;

  double wLo = set[iLo]
                 ? set[iLo]->GetVariation(0,set[iLo]->GetBin({q0,q3})) : 1.0;
  double wHi = set[iHi]
                 ? set[iHi]->GetVariation(0,set[iHi]->GetBin({q0,q3})) : 1.0;

  double t = (E - fEnergies[iLo]) / (fEnergies[iHi] - fEnergies[iLo]);
  double w = (1. - t) * wLo + t * wHi;
  return std::clamp(w, fWeightLimits.first, fWeightLimits.second);
}

//--------------------------------------------------------------------
// main event-response
event_unit_response_t
ValenciaMECq0q3InterpWeighting::GetEventResponse(genie::EventRecord const &ev)
{
  event_unit_response_t r;
  auto const &meta = GetSystMetaData();
  if (meta.empty()) return r;
  paramId_t pid = meta.front().systParamId;
  r.push_back({pid,{}});   // create the slot

  // apply only to CC νμ MEC on argon
  if (!ev.Summary()->ProcInfo().IsWeakCC() ||
    !ev.Summary()->ProcInfo().IsMEC()     ||
    (ev.Probe()->Pdg()!=genie::kPdgNuMu && ev.Probe()->Pdg()!=genie::kPdgAntiNuMu) ||
    (ev.Summary()->InitState().Tgt().Z()!=18) )
  {
    r.back().responses = {1.,1.,1.};
    return r;
  }

  double q0,q3; computeQ0Q3(ev,q0,q3);
  double E = ev.Probe()->P4()->E();
  double wCV = getInterpolatedWeight(E,q0,q3,classifyEvent(ev));

  double wUp = wCV;
  double wDn = (wCV>0) ? (2. - wCV) : 1.0;
  r.back().responses = {wDn,wCV,wUp};
  return r;
}
//--------------------------------------------------------------------
} // namespace nusyst
