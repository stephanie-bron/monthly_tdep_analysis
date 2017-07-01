#include "llh/public/CoordinateTransform.h"

#include "TMath.h"

#include "coord_interface/public/SLAHeader.h"


void LocalToEq(double zenDeg, double azDeg, double timeMJD,
	       double& raDeg, double & decDeg)
{
  double zenRad = zenDeg * TMath::DegToRad();
  double azRad  = azDeg  * TMath::DegToRad();

  raDeg = TMath::RadToDeg() * 
    SLACoordinateTransform::Local2RA( zenRad, azRad, timeMJD);
  decDeg = TMath::RadToDeg() * 
    SLACoordinateTransform::Local2Dec(zenRad, azRad, timeMJD);
}

void EqToLocal(double raDeg, double decDeg, double timeMJD,
	       double& zenDeg, double& azDeg)
{
  double raRad  = raDeg  * TMath::DegToRad();
  double decRad = decDeg * TMath::DegToRad();

  zenDeg = TMath::RadToDeg() * 
    SLACoordinateTransform::Equa2LocalZenith( raRad, decRad, timeMJD);
  azDeg  = TMath::RadToDeg() * 
    SLACoordinateTransform::Equa2LocalAzimuth(raRad, decRad, timeMJD);
}

void LocalToGalactic(double zenDeg, double azDeg, double timeMJD,
           double& galLonDeg, double & galLatDeg)
{
  double zenRad = zenDeg * TMath::DegToRad();
  double azRad  = azDeg  * TMath::DegToRad();

  galLonDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Local2GalacticLongitude(zenRad, azRad, timeMJD);
  galLatDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Local2GalacticLatitude(zenRad, azRad, timeMJD);
}

void GalacticToLocal(double galLonDeg, double galLatDeg, double timeMJD,
           double& zenDeg, double& azDeg)
{
  double galLonRad = galLonDeg * TMath::DegToRad();
  double galLatRad = galLatDeg * TMath::DegToRad();

  zenDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Galactic2LocalZenith(galLonRad, galLatRad, timeMJD);
  azDeg  = TMath::RadToDeg() *
    SLACoordinateTransform::Galactic2LocalAzimuth(galLonRad, galLatRad, timeMJD);
}

void EqToGalactic(double raDeg, double decDeg, double timeMJD,
           double& galLonDeg, double& galLatDeg)
{
  double raRad  = raDeg  * TMath::DegToRad();
  double decRad = decDeg * TMath::DegToRad();

  galLonDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Equatorial2GalacticLongitude( raRad, decRad, timeMJD);
  galLatDeg  = TMath::RadToDeg() *
    SLACoordinateTransform::Equatorial2GalacticLatitude(raRad, decRad, timeMJD);
}

void GalacticToEq(double galLonDeg, double galLatDeg, double timeMJD,
           double& raDeg, double & decDeg)
{
  double galLonRad = galLonDeg * TMath::DegToRad();
  double galLatRad = galLatDeg * TMath::DegToRad();

  raDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Galactic2RA( galLonRad, galLatRad, timeMJD);
  decDeg = TMath::RadToDeg() *
    SLACoordinateTransform::Galactic2Dec(galLonRad, galLatRad, timeMJD);
}





