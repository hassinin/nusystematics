#include "nusystematics/systproviders/QEInterference_tool.hh"
#include "nusystematics/utility/exceptions.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "Framework/GHEP/GHepParticle.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

QEInterference::QEInterference(ParameterSet const & ps) :
  IGENIESystProvider_tool(ps),
  QEIntfResponseCalculator(nullptr),
  ResponseParameterIdx(kParamUnhandled<size_t>) {}

QEInterference::~QEInterference() {}

bool QEInterference::SetupResponseCalculator(ParameterSet const & tool_options)
{
  // Silence GENIE
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");
  //genie::Messenger::Instance()->SetPriorityLevel("GHepUtils",
  //                                               log4cpp::Priority::FATAL);

  // Check the metadata makes sense
  if (!HasParam(GetSystMetaData(), "QEInterference")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("QEInterference");
  }

  // Get manifests for options
  if (!tool_options.has_key("QEInterference_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: QEInterference parameter exists in the "
           "SystMetaData, but "
           "no QEInterference_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "QEInterference::BuildSystMetaData, but wasn't, this is a "
           "bug, "
           "please report to the maintainer.";
  }

  fhicl::ParameterSet const &templateManifest =
      tool_options.get<fhicl::ParameterSet>(
          "QEInterference_input_manifest");

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "QEInterference");

  // Initialise the calculator
  QEIntfResponseCalculator = std::make_unique<QEInterferenceResponseCalculator>(templateManifest);

  return true;
}

SystMetaData QEInterference::BuildSystMetaData(ParameterSet const & cfg,
						     paramId_t id) 
{
  SystMetaData smd;

  SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "QEInterference",
                                                 phdr, id)) {
    phdr.systParamId = id++;
    smd.push_back(phdr);
  }

  ParameterSet templateManifest =
      cfg.get<ParameterSet>("QEInterference_input_manifest");

  if (!cfg.has_key("QEInterference_input_manifest") ||
      !cfg.is_key_to_table("QEInterference_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for "
           "QEInterference, expected to find a FHiCL table keyed by "
           "QEInterference_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }

  tool_options.put("QEInterference_input_manifest", templateManifest);

  return smd;
}

std::string QEInterference::AsString() { return "QEInterference"; }

event_unit_response_t QEInterference::GetEventResponse(genie::EventRecord const & ev) 
{ 
  // Anything but QE is not handled and should be default
  // Use GetDefaultEventResponse() to return an auto-1.-filled vector

  if ( !ev.Summary()->ProcInfo().IsQuasiElastic() ||
       (!ev.Summary()->ProcInfo().IsWeakCC() &&
	!ev.Summary()->ProcInfo().IsWeakNC()) ) {
    return this->GetDefaultEventResponse();
  }

  // Load up the Enu, Q0, and Q3. Also the PDG code!
  genie::GHepParticle * neutrino = ev.Probe();
  genie::GHepParticle * FSLep    = ev.FinalStatePrimaryLepton();

  if( ! neutrino || ! FSLep ) {
    throw incorrectly_generated()
      << "[ERROR]: Failed to find neutrino and final-state primary lepton in event: "
      << ev.Summary()->AsString();
  }

  TLorentzVector p4_nu  = *( neutrino->GetP4() );
  double Enu = p4_nu.E();
  int    pdg = neutrino->Pdg();
  int current = static_cast<int>( ev.Summary()->ProcInfo().IsWeakNC() ); // 0 if CC, 1 if NC
  
  TLorentzVector p4_lep = *( FSLep->GetP4() );
  
  // Q0 == Enu - Elep
  // Q3 == ||\vec{p3nu} - \vec{p3lep}||
  TLorentzVector delta_p4 = p4_nu - p4_lep;
  double Q0 = delta_p4.E();
  double Q3 = delta_p4.P();

  // Make the output
  event_unit_response_t resp;
  SystParamHeader const & hdr = GetSystMetaData()[ResponseParameterIdx];

  resp.push_back( {hdr.systParamId, {}} );
  for (double var : hdr.paramVariations) {
    double this_reweight = (1-var) + var * QEIntfResponseCalculator->GetWeight( pdg, current,
										Enu, Q0, Q3 );
    resp.back().responses.push_back( this_reweight );
  }

  return resp;
}
