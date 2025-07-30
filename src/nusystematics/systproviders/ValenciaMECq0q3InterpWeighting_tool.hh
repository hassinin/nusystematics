
#pragma once
#include "nusystematics/utility/exceptions.hh"
#include "fhiclcpp/ParameterSet.h"
#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/responsecalculators/ValenciaMECq0q3ResponseCalc.hh"
#include "systematicstools/interface/EventResponse_product.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/interface/ISystProviderTool.hh"
#include "EventRecord.h"
#include "TH2.h"
#include "TFile.h"
#include <map>
#include <vector>
#include <memory>

namespace nusyst {

enum class Topology { np, nn, unknown };

class ValenciaMECq0q3InterpWeighting : public IGENIESystProvider_tool {
public:
    ValenciaMECq0q3InterpWeighting(fhicl::ParameterSet const &ps);
    systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &cfg, systtools::paramId_t firstId);
    bool SetupResponseCalculator(fhicl::ParameterSet const &tool_opts);
    fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }
    Topology classifyEvent(genie::EventRecord const &ev) const;
    void computeQ0Q3(genie::EventRecord const &ev, double &q0, double &q3) const;
    double getInterpolatedWeight(double E, double q0, double q3, Topology topo) const;
    systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &ev);

private:
    std::vector<double> fEnergies;
    std::map<Topology, std::vector<std::unique_ptr<ValenciaMECq0q3ResponseCalc>>> fCalcs;
    std::pair<double, double> fWeightLimits{0.1, 5.0};
    fhicl::ParameterSet tool_options;
};

} // namespace nusyst
