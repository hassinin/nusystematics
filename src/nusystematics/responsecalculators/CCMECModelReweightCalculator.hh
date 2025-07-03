#ifndef nusystematics_RESPONSE_CALCULATORS_CCMECModelReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCMECModelReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

namespace nusyst {

  class CCMECModelReweightCalculator {
  public:
    // Load input histograms, given the filename
    void LoadInputHistograms(const std::string &file);

    // Calculate the reweight, given relevant kienatmics and the sigma of the variation
    double GetReweight(double Enu, double q0, double q3, int nn_np, double sigma);

    // Return the name of the calculator
    std::string GetCalculatorName() const { return "CCMECModelReweightCalculator"; }

    // class constructor
    CCMECModelReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest.get<std::string>("input_file"));
    }

  protected:
    // reweight histograms for NN interactions
    std::vector<std::unique_ptr<TH2D>> fReweightNN;

    // reweight histograms for NP interactions
    std::vector<std::unique_ptr<TH2D>> fReweightNP;
  };

} // end namespace nusyst

inline void nusyst::CCMECModelReweightCalculator::LoadInputHistograms(const std::string &file) {
  // TODO: implement
}

inline double nusyst::CCMECModelReweightCalculator::GetReweight(double Enu, double q0, double q3, int nn_np, double sigma) {
  // TODO: implement
  return 1.;
}

#endif
