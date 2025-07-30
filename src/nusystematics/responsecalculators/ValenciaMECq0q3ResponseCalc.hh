#pragma once
//
// Simple 2-D (q0,q3) response calculator for Valencia MEC.
//
// It inherits from TemplateResponseCalculatorBase<2,false> (2-D, NOT
// vertical) so that NuSystematics handle() helpers will work.
//
#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

#include "TH2.h"

namespace nusyst {

class ValenciaMECq0q3ResponseCalc : public TemplateResponseCalculatorBase<2,false> {
public:
  using base_t = TemplateResponseCalculatorBase<2,false>;

  ValenciaMECq0q3ResponseCalc() = default;
  /// Construct from an *already-cloned* TH2 histogram
  explicit ValenciaMECq0q3ResponseCalc(TH2 *hist) {
    // Store the histogram in BinnedResponses with key 0
    this->BinnedResponses[0] = std::unique_ptr<TH2>(hist);
  }

  std::string GetCalculatorName() const override;

  /** Return the (q0,q3) bin index with edge clamping. */
  bin_it_t GetBin(std::array<double,2> const &q0q3) const override {
    TH2 *h = this->BinnedResponses.begin()->second.get();
    int x = h->GetXaxis()->FindFixBin(q0q3[1]);          // q3 on X
    int y = h->GetYaxis()->FindFixBin(q0q3[0]);          // q0 on Y
    x = std::clamp(x, 1, h->GetNbinsX());
    y = std::clamp(y, 1, h->GetNbinsY());
    return h->GetBin(x,y);
  }
};

  // ...existing code...

} // namespace nusyst
