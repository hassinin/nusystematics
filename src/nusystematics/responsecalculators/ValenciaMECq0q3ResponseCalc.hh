/*******************************************************************************
 * 1.  src/ValenciaMECq0q3ResponseCalc.hh
 ******************************************************************************/
#ifndef NUSYST_VALENCIA_MEC_Q0Q3_RESPONSECALC_HH
#define NUSYST_VALENCIA_MEC_Q0Q3_RESPONSECALC_HH

#include <memory>
#include <TH2D.h>

namespace nusyst {

/// Lightweight helper that evaluates a single 2‑D weight histogram.
class ValenciaMECq0q3ResponseCalc {
public:
  /// Constructor takes ownership of the histogram (cloned internally).
  ValenciaMECq0q3ResponseCalc(TH2D* h, double w_min = 0.0, double w_max = 5.0);

  /// Central weight with **bilinear interpolation** inside the map.
  double GetCentralWeight(double q0, double q3) const;

  /// Up / Down side variation (ivar = ±1) – symmetric envelope 2 − w_CV.
  double GetVariation(int ivar, double q0, double q3) const;

private:
  std::unique_ptr<TH2D> fHist;   ///< owned, thread‑safe clone of the histogram
  double fWmin, fWmax;           ///< clamp limits
};

} // namespace nusyst

#endif // NUSYST_VALENCIA_MEC_Q0Q3_RESPONSECALC_HH