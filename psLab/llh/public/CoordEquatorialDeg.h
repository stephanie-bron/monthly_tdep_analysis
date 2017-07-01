#ifndef LLH_COORDEQUATORIALDEG_H_
#define LLH_COORDEQUATORIALDEG_H_

#include "TMath.h"
#include "TVector3.h"

#include "llh/public/CoordClasses.h"


class Spherical : public Coord {
  TVector3 v_;
public:
  Spherical () : v_(0.,0., 1.) {}

  void SetCoordsPolarRad (double thetaRad, double phiRad) 
  { v_.SetMagThetaPhi(1.,thetaRad,phiRad); }

  double GetThetaRad() const {return v_.Theta(); }
  double GetPhiRad() const {return v_.Phi(); }

  void SetCoordsPolarDeg (double thetaDeg, double phiDeg) 
  { SetCoordsPolarRad(thetaDeg*TMath::DegToRad(), 
		      phiDeg*TMath::DegToRad() ); }
  double GetThetaDeg() const {return v_.Theta()*TMath::RadToDeg(); }
  double GetPhiDeg() const {return v_.Phi()*TMath::RadToDeg(); }


  double DistanceTo (const Coord& coord2) const
  {
    const Spherical* s2 = 
      dynamic_cast<const Spherical*>(&coord2);
    assert(s2);
    return v_.Angle( s2->v_);
  }


  void RotateX(double angleRad) {
    v_.RotateX(angleRad);
  }
  void RotateY(double angleRad) {
    v_.RotateY(angleRad);
  }
  void RotateZ(double angleRad) {
    v_.RotateZ(angleRad);
  }
};



class EquatorialDeg : public Coord {
  double raDeg_;
  double decDeg_;
  mutable bool s_IsSet_;  // save time calculating s_ until it is needed
  mutable Spherical s_;

  void Set_s() const {
      s_.SetCoordsPolarDeg(90.-decDeg_, raDeg_); 
      s_IsSet_ = true;  // so won't have to do this again if coords unchanged
  }

public:
  EquatorialDeg () : raDeg_(0.), decDeg_(90.), s_IsSet_(false) {}

  EquatorialDeg(double raDeg, double decDeg) { SetCoords(raDeg, decDeg);}

  void SetCoords (double raDeg, double decDeg) {
    raDeg_ = raDeg;
    decDeg_ = decDeg;
    s_IsSet_ = false;  // coords have changed, s_ not recalculated yet
  }

  void SetRaDeg(double raDeg) { 
    raDeg_ = raDeg;
    s_IsSet_ = false;  // coords have changed, s_ not recalculated yet
  }

  void SetDecDeg(double decDeg) {
    decDeg_ = decDeg;
    s_IsSet_ = false;  // coords have changed, s_ not recalculated yet
  }


  double GetRa() const {return raDeg_; }
  double GetDec() const {return decDeg_; }

  double DistanceTo (const Coord& coord2) const
  {
    const EquatorialDeg* e2 = 
      dynamic_cast<const EquatorialDeg*>(&coord2);
    assert(e2);

    if (!     s_IsSet_) {     Set_s(); }
    if (! e2->s_IsSet_) { e2->Set_s(); }

    return s_.DistanceTo(e2->s_)*TMath::RadToDeg();
  }

};


#endif // LLH_COORDEQUATORIALDEG_H_
