#include "llh/public/CoordClasses.h"

#include "TMath.h"


double CircularGaussUnc(double r, double sigma) {
  return exp(-r*r/(sigma*sigma*2)) / (2.*TMath::Pi()*sigma*sigma);
}
