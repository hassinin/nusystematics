#ifndef nusystematics_SYSTPROVIDERS_CCMECModelReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_CCMECModelReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/CCMECModelReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class CCMECModelReweight: public nusyst::IGENIESystProvider_tool {

  // Constructor
  explicit CCMECModelReweight(const fhicl::ParameterSet &ps):
    IGENIESystProvider_tool(ps),
    fCalculator(ps) {}

  // Calculate weight
  systtools::event_unit_response_t GetEventResponse(const genie::EventRecord &ev);

protected:
  nusyst::CCMECModelReweightCalculator fCalculator;
};
#endif
