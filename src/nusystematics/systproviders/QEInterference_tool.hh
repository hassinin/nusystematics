#pragma once

// Interfaces
#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "nusystematics/responsecalculators/QEInterferenceResponseCalculator.hh"

#include "Framework/Messenger/Messenger.h"

#include "TFile.h"
#include "TH3D.h"

#include <string>
#include <sstream>

/*
 * \class QEInterference
 *
 * \brief A weighting tool that applies the ratio of ACHILLES-derived cross sections:
 *        --> Numerator: Xsec with the 1p1h + 2p2h interference term
 *        switched ON for 1-particle knockout events
 *        --> Denominator: Xsec with the interference term switched OFF.
 *
 * \reference arXiv: 2312.12545 
 *
 * \created September 22, 2025
 *
 * \authors John Plows <kplows@liverpool.ac.uk>
 *          Gray Putnam <putnam@fnal.gov>
 */

class QEInterference : public nusyst::IGENIESystProvider_tool {

  std::unique_ptr<nusyst::QEInterferenceResponseCalculator> QEIntfResponseCalculator;

public:
  explicit QEInterference(fhicl::ParameterSet const & ps);
  ~QEInterference();
  
  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const & cfg,
					    systtools::paramId_t id);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);
  
  std::string AsString();

private:
  
  fhicl::ParameterSet tool_options;

  // Vector to hold the name of the descriptors.
  // Will be loaded into the appropriate members in BuildSystMetaData()
  std::vector<std::string> descriptors = {
    "QEIntf_dial", ///< Tweak dial t controlling linear response: f(t; W) = 1-t + t*W
    "ACHILLES_strength" ///< Linear scaling to weight. g(a; W) = a*W
  };

  //size_t ResponseParameterIdx;
  std::vector<size_t> ResponseParameterIndices;

  // Members containing the tweak dials, and tweak strength
  //std::vector<std::vector<double>> variations;

  int nu_pdg;
  double Enu, Q0, Q3;

  int verbosity_level;

}; // class QEInterference
