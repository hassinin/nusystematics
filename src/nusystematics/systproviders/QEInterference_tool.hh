#pragma once

// Interfaces
#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "nusystematics/responsecalculators/QEInterferenceResponseCalculator.hh"

#include "Framework/Messenger/Messenger.h"

#include "TFile.h"
#include "TH3D.h"

#include <string>

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

  size_t ResponseParameterIdx;

  int nu_pdg;
  double Enu, Q0, Q3;

}; // class QEInterference
