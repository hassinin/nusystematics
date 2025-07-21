#ifndef nusystematics_SYSTPROVIDERS_CCMECModelReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_CCMECModelReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/CCMECModelReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <string>

// Forward declaration of genie::EventRecord
namespace genie {
class EventRecord;
}

namespace nusyst {
class CCMECModelReweight : public IGENIESystProvider_tool {

public:
  // Constructor from a FHiCL parameter set
  explicit CCMECModelReweight(fhicl::ParameterSet const &ps);

  // The primary method for calculating the event weight
  systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &ev) override;

protected:
  // The calculator object that loads histograms and performs the reweighting logic
  CCMECModelReweightCalculator fCalculator;
};
} // namespace nusyst

#endif
