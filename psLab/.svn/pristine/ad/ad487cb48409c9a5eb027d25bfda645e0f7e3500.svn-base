#ifndef LLH_COORDINATE_TRANSFORM_H_
#define LLH_COORDINATE_TRANSFORM_H_

// as fraction of Julian (Solar) Day
#define LENGTH_OF_SIDEREAL_DAY_FRAC  0.997270
#define LENGTH_OF_SIDEREAL_DAY_SECONDS 86164.128

void LocalToEq(double zenDeg, double azDeg, double timeMJD,
	       double& raDeg, double & decDeg);

void EqToLocal(double raDeg, double decDeg, double timeMJD,
	       double& zenDeg, double& azDeg);

void LocalToGalactic(double zenDeg, double azDeg, double timeMJD,
           double& galLonDeg, double & galLatDeg);

void GalacticToLocal(double galLonDeg, double galLatDeg, double timeMJD,
           double& zenDeg, double& azDeg);

void EqToGalactic(double raDeg, double decDeg, double timeMJD,
           double& galLonDeg, double& galLatDeg);

void GalacticToEq(double galLonDeg, double galLatDeg, double timeMJD,
           double& raDeg, double & decDeg);


double LengthOfSiderealDayInSeconds() { return LENGTH_OF_SIDEREAL_DAY_SECONDS; }

#endif // LLH_COORDINATE_TRANSFORM_H_
