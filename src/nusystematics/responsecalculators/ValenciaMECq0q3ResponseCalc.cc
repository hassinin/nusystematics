/*******************************************************************************
 * 2.  src/ValenciaMECq0q3ResponseCalc.cc
 ******************************************************************************/
#include "ValenciaMECq0q3ResponseCalc.hh"
#include <algorithm>
#include <stdexcept>

using namespace nusyst;

// ---------------------------------------------------------------------------
ValenciaMECq0q3ResponseCalc::ValenciaMECq0q3ResponseCalc(TH2D* h,
                                                         double wmin,
                                                         double wmax)
  : fHist(h ? static_cast<TH2D*>(h->Clone()) : nullptr),
    fWmin(wmin), fWmax(wmax)
{
  if (!fHist)
    throw std::runtime_error("Null histogram passed to ResponseCalc");
  fHist->SetDirectory(nullptr);          // detach from any ROOT file
  // Note: kCanRebin is not available in this ROOT version, binning is fixed by design
}

// ---------------------------------------------------------------------------
inline double clamp(double x, double a, double b) {
  return std::max(a, std::min(b, x));
}

// ---------------------------------------------------------------------------
double ValenciaMECq0q3ResponseCalc::GetCentralWeight(double q0, double q3) const
{
  // Histogram axis order = (x = q3, y = q0)
  const double xMin = fHist->GetXaxis()->GetXmin();
  const double xMax = fHist->GetXaxis()->GetXmax();
  const double yMin = fHist->GetYaxis()->GetXmin();
  const double yMax = fHist->GetYaxis()->GetXmax();

  // If kinematics fall outside trained / provided map domain, return neutral weight
  if (q3 < xMin || q3 > xMax || q0 < yMin || q0 > yMax) {
    return 1.0; // do not extrapolate beyond Valencia map coverage
  }

  const double w = fHist->Interpolate(q3, q0); // bilinear inside domain
  return clamp(w, fWmin, fWmax);
}

// ---------------------------------------------------------------------------
double ValenciaMECq0q3ResponseCalc::GetVariation(int ivar,
                                                 double q0, double q3) const
{
  if (ivar == 0) return GetCentralWeight(q0, q3);
  const double w_cv = GetCentralWeight(q0, q3);
  return clamp(2.0 - w_cv, fWmin, fWmax);  // symmetric envelope
}