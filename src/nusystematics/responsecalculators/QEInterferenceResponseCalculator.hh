#pragma once

/*
 * \class QEInterferenceResponseCalculator
 *
 * \brief Looks up from a histogram in (Enu, Q0, Q3) space the bin for an event, and applies it.
 *        Based on CCQERPAReweightCalculator.hh for boilerplate code
 *
 * \created September 22, 2025
 *
 * \authors John Plows <kplows@liverpool.ac.uk>
 *          Gray Putnam <putnam@fnal.gov>
 */

#include <algorithm>
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"
#include "fhiclcpp/ParameterSet.h"
#include "TH1D.h"
#include "TH3D.h"

NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_INPUTHIST);
NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_FILEPATH);
NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_BIN);

namespace nusyst {

  class QEInterferenceResponseCalculator {

  public: 
    QEInterferenceResponseCalculator(fhicl::ParameterSet const & InputParams) {
      LoadInputHistograms(InputParams);
    }

    ~QEInterferenceResponseCalculator() {}

    void LoadInputHistograms(fhicl::ParameterSet const & InputParams);

    double GetWeight(int pdg, double Enu_GeV, double Q0_GeV, double Q3_GeV);

    std::string GetCalculatorName() const { return "QEInterferenceResponseCalculator"; }

  private:

    // RETHERE: different bins per flavour?
    std::vector<double> Enu_bins, Q0_bins, Q3_bins;
    TH3D numu_ratio_histogram;
    
  }; // class QEinterferenceResponseCalculator

  inline double QEInterferenceResponseCalculator::GetWeight( int pdg, double Enu_GeV,
							     double Q0_GeV, double Q3_GeV ) {
    // If outside the kinematic range, return 1.0
    if( Enu_GeV < *(Enu_bins.begin()) || Enu_GeV > *(Enu_bins.end()-1) ||
	Q0_GeV < *(Q0_bins.begin()) || Q0_GeV > *(Q0_bins.end()-1) ||
	Q3_GeV < *(Q3_bins.begin()) || Q3_GeV > *(Q3_bins.end()-1) ) return 1.0;

    // Seek out the bin and return it

    std::vector<double>::iterator it_enu = Enu_bins.begin();
    std::vector<double>::iterator it_q0  = Q0_bins.begin();
    std::vector<double>::iterator it_q3  = Q3_bins.begin();
    
    while( it_enu < Enu_bins.end() && Enu_GeV >= *it_enu ){ ++it_enu; } --it_enu;
    while( it_q0 < Q0_bins.end() && Q0_GeV >= (*it_q0) ){ ++it_q0; } --it_q0;
    while( it_q3 < Q3_bins.end() && Q3_GeV >= (*it_q3) ){ ++it_q3; } --it_q3;

    int bin_enu = it_enu - Enu_bins.begin();
    int bin_q0  = it_q0  -  Q0_bins.begin();
    int bin_q3  = it_q3  -  Q3_bins.begin();
    switch(pdg) {
    case 14:
      return numu_ratio_histogram.GetBinContent( bin_enu, bin_q0, bin_q3 );
    default: // RETHERE: to add
      return 1.0;
    } // switch pdg
  } // GetWeight()

  inline void QEInterferenceResponseCalculator::LoadInputHistograms(fhicl::ParameterSet const & ps)
  {
    const std::string & default_input_file = ps.get<std::string>("input_file", "");
    const std::vector<std::string> known_flavours = { "numu", "numubar", "nue", "nuebar" };
    
    // Obtain for each flavour the name of the input histogram, and load it in.
    for( fhicl::ParameterSet const & val_config :
	   ps.get<std::vector<fhicl::ParameterSet>>("inputs") ) {
      std::string flavour_name = val_config.get<std::string>("name"); // for flavours. nu(mu,e)(bar)
      std::string input_hist   = val_config.get<std::string>("input_hist"); // name of the hist in file
      std::string input_file   = val_config.get<std::string>("input_file", default_input_file);

      if( std::find( known_flavours.begin(), known_flavours.end(), flavour_name.c_str() ) ==
	  known_flavours.end() ) {
	throw invalid_QEIntf_INPUTHIST() << "[ERROR]: Unknown input flavour " << flavour_name.c_str()
					 << " ; check inputs:name";
      }

      // If input_file is not given as an absolute path seatch ${NUSYSTEMATICS_FQ_DIR}/data/
      if( input_file.find("/") != 0 ) {
	if( std::getenv("nusystematics_ROOT") == "" ) {
	  throw invalid_QEIntf_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set!\n"
					  << "Given path: " << input_file << "\n"
					  << "Expect absolute path (starts with '/'), or file relative to ${nusystematics_ROOT}/data/";
	}
	input_file = std::string(std::getenv("nusystematics_ROOT")) + "/data/" + input_file;
      } // absolute path not given

      if( flavour_name == "numu" ){
	numu_ratio_histogram = *(GetHistogram<TH3D>( input_file, input_hist ));
      }
      TH1D hEnu_bins = *(GetHistogram<TH1D>( input_file, "hEnuBins" ));
      TH1D hQ0_bins  = *(GetHistogram<TH1D>( input_file, "hQ0Bins" ));
      TH1D hQ3_bins  = *(GetHistogram<TH1D>( input_file, "hQ3Bins" ));

      for( Int_t ib = 1; ib <= hEnu_bins.GetNbinsX(); ib++ )
	Enu_bins.emplace_back( hEnu_bins.GetBinLowEdge(ib) );
      for( Int_t ib = 1; ib <= hQ0_bins.GetNbinsX(); ib++ ) 
	Q0_bins.emplace_back( hQ0_bins.GetBinLowEdge(ib) );
      for( Int_t ib = 1; ib <= hQ3_bins.GetNbinsX(); ib++ ) 
	Q3_bins.emplace_back( hQ3_bins.GetBinLowEdge(ib) );

      // Add the high edges
      Enu_bins.emplace_back( hEnu_bins.GetBinLowEdge(hEnu_bins.GetNbinsX()) + 
			     hEnu_bins.GetBinWidth(hEnu_bins.GetNbinsX()) );
      Q0_bins.emplace_back( hQ0_bins.GetBinLowEdge(hQ0_bins.GetNbinsX()) + 
			    hQ0_bins.GetBinWidth(hQ0_bins.GetNbinsX()) );
      Q3_bins.emplace_back( hQ3_bins.GetBinLowEdge(hQ3_bins.GetNbinsX()) + 
			    hQ3_bins.GetBinWidth(hQ3_bins.GetNbinsX()) );

    } // for each flavour in inputs
  } // LoadInputHistograms

} // namespace nusyst
