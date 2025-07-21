#include "nusystematics/systproviders/CCMECModelReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/utility/GENIEUtils.hh" // For nusyst::pdg
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/ProcInfo.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"

#include "TLorentzVector.h"

#include "cetlib/PluginTypeDeducer.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

// Constructor implementation
CCMECModelReweight::CCMECModelReweight(fhicl::ParameterSet const &ps)
    : IGENIESystProvider_tool(ps), fCalculator(ps) {}

// Main implementation of the event reweighting
event_unit_response_t
CCMECModelReweight::GetEventResponse(genie::EventRecord const &ev) {

  // --- Pre-selection for performance ---
  if (!ev.Summary()->ProcInfo().IsMEC() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return GetDefaultEventResponse();
  }

  if (ev.Summary()->InitState().Tgt().A() != 40) {
    return GetDefaultEventResponse();
  }

  // Explicitly check for muon neutrino (PDG code 14)
  if (ev.Summary()->InitState().ProbePdg() != genie::kPdgNuMu) {
    return GetDefaultEventResponse();
  }

  // --- Extract Kinematics ---
  const genie::Interaction *interaction = ev.Summary();
  double Enu = interaction->InitState().ProbeE(genie::kRfLab);
  const TLorentzVector &k_nu = *interaction->InitState().GetProbeP4();
  const TLorentzVector &k_lep = *interaction->FinalState().GetLeptonP4();
  TLorentzVector q = k_nu - k_lep;
  double q0 = q.E();
  double q3 = q.Vect().Mag();

  // --- Determine Initial State by Inspecting the FINAL STATE Nucleons ---
  // This logic is specific to numu CC interactions (v_mu + n -> mu- + p)
  // and matches the binning of my ROOT file.
  int n_final_protons = 0;
  int n_final_neutrons = 0;
  for (const genie::GHepParticle &p : ev) {
    if (p.Status() == genie::kIStStableFinalState && p.IsNucleon()) {
      if (p.Pdg() == genie::kPdgProton) {
        n_final_protons++;
      } else if (p.Pdg() == genie::kPdgNeutron) {
        n_final_neutrons++;
      }
    }
  }

  int nn_np_id = -1;
  // For a numu CC MEC event:
  // An initial 'nn' pair results in a 'pn' final state (v + nn -> mu- + p + n)
  if (n_final_protons == 1 && n_final_neutrons == 1) {
    nn_np_id = 1; // This event came from an initial nn-pair
  }
  // An initial 'np' pair results in a 'pp' final state (v + np -> mu- + p + p)
  else if (n_final_protons == 2 && n_final_neutrons == 0) {
    nn_np_id = 0; // This event came from an initial np-pair
  }
  else {
    // If the final state doesn't have exactly two nucleons, we can't
    // reliably identify it for this reweighting. Return a neutral weight.
    return GetDefaultEventResponse();
  }


  // --- Calculate Weight ---
  double sigma = GetSystParams(0)[0];
  double weight = fCalculator.GetReweight(Enu, q0, q3, nn_np_id, sigma);

  // --- Package and Return Response ---
  event_unit_response_t resp;
  resp.push_back({GetParIDs(0)[0], weight});
  return resp;
}

// Register this tool with the ART/LArSoft framework
DEFINE_ART_PLUGIN(CCMECModelReweight)
