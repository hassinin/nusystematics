#ifndef nusystematics_RESPONSE_CALCULATORS_CCMECModelReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCMECModelReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TFile.h"
#include "TH2.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace nusyst {

class CCMECModelReweightCalculator {
public:
  void LoadInputHistograms(const std::string &file);
  double GetReweight(double Enu, double q0, double q3, int nn_np, double sigma);

  std::string GetCalculatorName() const {
    return "CCMECModelReweightCalculator";
  }

  explicit CCMECModelReweightCalculator(fhicl::ParameterSet const &ps) {
    LoadInputHistograms(ps.get<std::string>("input_file"));
  }

protected:
  // reweight histograms for NN interactions
  std::vector<std::unique_ptr<TH2D>> fReweightNN;
  // reweight histograms for NP interactions
  std::vector<std::unique_ptr<TH2D>> fReweightNP;
  // Neutrino energy points corresponding to the histograms
  std::vector<double> fEnergies;
};

} // namespace nusyst

inline void nusyst::CCMECModelReweightCalculator::LoadInputHistograms(
    const std::string &file) {

  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) {
    throw systtools::exceptions::CouldNotOpenFile(file);
  }

  const double E_min = 0.4;
  const double E_max = 2.5;
  const double E_step = 0.1;
  const int n_steps = static_cast<int>((E_max - E_min) / E_step + 1.5);

  for (int i = 0; i < n_steps; ++i) {
    double energy = E_min + i * E_step;
    fEnergies.push_back(energy);

    // --- THIS IS THE MODIFIED SECTION ---

    // Construct names for both np and nn histograms
    TString hist_name_np = TString::Format("h_weights_map_np_%.1fGeV", energy);
    hist_name_np.ReplaceAll(".", "p");

    TString hist_name_nn = TString::Format("h_weights_map_nn_%.1fGeV", energy);
    hist_name_nn.ReplaceAll(".", "p");

    // Load the NP histogram
    TH2D *h_np = systtools::root::LoadObject<TH2D>(&f, hist_name_np.Data(), true);
    h_np->SetDirectory(nullptr); // Detach from file
    fReweightNP.push_back(std::unique_ptr<TH2D>(h_np));

    // Load the NN histogram
    TH2D *h_nn = systtools::root::LoadObject<TH2D>(&f, hist_name_nn.Data(), true);
    h_nn->SetDirectory(nullptr); // Detach from file
    fReweightNN.push_back(std::unique_ptr<TH2D>(h_nn));
  }
  f.Close();
}

inline double
nusyst::CCMECModelReweightCalculator::GetReweight(double Enu, double q0,
                                                  double q3, int nn_np,
                                                  double sigma) {

  // This logic selects between the different nn and np histograms
  const auto &hists = (nn_np == 1) ? fReweightNN : fReweightNP;

  if (Enu < fEnergies.front() || Enu > fEnergies.back()) {
    return 1.0;
  }

  auto it =
      std::upper_bound(fEnergies.begin(), fEnergies.end(), Enu);
  auto high_E_it = it;
  auto low_E_it = std::prev(it);

  size_t high_E_idx = std::distance(fEnergies.begin(), high_E_it);
  size_t low_E_idx = std::distance(fEnergies.begin(), low_E_it);

  const auto &h_low = hists.at(low_E_idx);
  const auto &h_high = hists.at(high_E_idx);
  
  int bin_low = h_low->FindBin(q3, q0);
  if (h_low->IsBinUnderflow(bin_low) || h_low->IsBinOverflow(bin_low)) {
    return 1.0;
  }
  double w_low = h_low->GetBinContent(bin_low);
  
  int bin_high = h_high->FindBin(q3, q0);
  if (h_high->IsBinUnderflow(bin_high) || h_high->IsBinOverflow(bin_high)) {
    return 1.0;
  }
  double w_high = h_high->GetBinContent(bin_high);
  
  double E_low = *low_E_it;
  double E_high = *high_E_it;
  double interpolated_weight = w_low + (w_high - w_low) * (Enu - E_low) / (E_high - E_low);
  
  if (interpolated_weight == 0) {
    return 1.0;
  }
  
  return 1.0 + sigma * (interpolated_weight - 1.0);
}

#endif