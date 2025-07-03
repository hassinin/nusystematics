#include "nusystematics/systproviders/CCMECModelReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

event_unit_response_t CCMECModelReweight::GetEventResponse(genie::EventRecord const &ev) {

  // Return default for non-MEC
  if (!ev.Summary()->ProcInfo().IsMEC() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return GetDefaultEventResponse();
  }


  // TODO: implement, and replace default
  return GetDefaultEventResponse();
}
