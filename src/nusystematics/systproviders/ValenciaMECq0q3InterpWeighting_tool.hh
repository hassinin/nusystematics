/*******************************************************************************
 *   ValenciaMECq0q3InterpWeighting_tool.hh
 ******************************************************************************/
#ifndef NUSYST_VALENCIA_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH
#define NUSYST_VALENCIA_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/responsecalculators/ValenciaMECq0q3ResponseCalc.hh"
#include <unordered_map>

namespace nusyst {

class ValenciaMECq0q3InterpWeighting : public IGENIESystProvider_tool {
public:
  explicit ValenciaMECq0q3InterpWeighting(const fhicl::ParameterSet& pset);
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);
  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const& ev) override;

private:
  fhicl::ParameterSet tool_options;
  enum class Topo { np = 0, nn = 1, unknown = 2 };
  static Topo  ClassifyEvent(genie::EventRecord const&);
  static void  ComputeQ0Q3(genie::EventRecord const&, double& q0, double& q3,
                           double& Enu);

  std::vector<double> fEgrid; ///< GeV
  double fWmin{0.0};
  double fWmax{5.0};
  std::unordered_map<Topo,
      std::vector<std::unique_ptr<ValenciaMECq0q3ResponseCalc>>> fCalcs;
};

} // namespace nusyst

#endif // NUSYST_VALENCIA_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH
