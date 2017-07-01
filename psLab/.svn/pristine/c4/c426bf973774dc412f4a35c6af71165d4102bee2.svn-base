#ifndef LLH_NEWLLHSTACK_H_
#define LLH_NEWLLHSTACK_H_

#include "TFitterMinuit.h"
#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

#include "llh/public/NewLlhEnergy.h"
#include "llh/public/CoordEquatorialDeg.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;
class NewLlhExtFCN;


class NewLlhExt : public NewLlhEnergy {

 private:
  // New for stacking/extended sources:
  NewLlhExtFCN* fcn_;
  double srcSigma_;
  vector<double> spaceRatioVects_;

  void OptimizeEventSelection();

 public:


  friend class NewLlhExtFCN;

  NewLlhExt();
  virtual ~NewLlhExt() { 
    if (minuit_) { 
      delete minuit_;
      minuit_=NULL;
    }
  }

  void PrepareAnalysis();

  // These are exact copies from NewLlhEnergy, but are required because fcn_ is of
  //  a different type (overloaded).  
  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double nSrc, double gamma);

  virtual double EvaluateLlh(double* parArray) {
    return EvaluateLlh(parArray[0], parArray[1]);
  }

  // New for Extended:
  /* @brief Copies external vector of source sigmas to interal member
  */
  void SetSourceSigma(double sourceSigma);

};



class NewLlhExtFCN : public ROOT::Minuit2::FCNBase {
 private:
  NewLlhExt* ptr;
 public:
  // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
  // This is what gets minimized: you have to define your likelihood
  virtual double operator() (const vector<double>& par) const;

  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<NewLlhExt*>(fn); 
  }
};

#endif // LLH_NEWLLHSTACK_H_
