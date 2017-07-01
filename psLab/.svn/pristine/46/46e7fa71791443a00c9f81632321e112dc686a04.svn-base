#ifndef LLH_NEWLLHTH2F_H_
#define LLH_NEWLLHTH2F_H_

#include "TFitterMinuit.h"
#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

#include "llh/public/NewLlhEnergy.h"
#include "llh/public/CoordEquatorialDeg.h"

#include "TH2F.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;
class NewLlhTH2FFCN;


class NewLlhTH2F : public NewLlhEnergy {

 private:
  // JD: We can still just use NewLlhEnergyFCN
  //NewLlhTH2FFCN* fcn_;


  void OptimizeEventSelection();


 public:
  // TODO: JD: public for debugging:
  TH2F hPdf_[30];

  // These are the reverse conversions as I3TH2FGenerator, to match input histos
  // Can converst Event coords in Equatorial to match input hist in gal coords
  bool convEqToGal_;
  // Event coords in deg can be converted to match input hist in radians
  bool convDegToRad_;
  // Event coords in Ra=Phi and Dec=Theta-90 can be shifted to match histos 
  //  in Theta and Phi
  bool shiftDecToTheta_;
  // If histo has galactic center at center of plot (I.e. binning range -180 to 
  //  180deg instead of 0 to 360deg), you have to shift here.  A corresponding
  //  shift in I3TH2FGenerator is not actually needed here.
  bool shiftGalLonBinning_;

  //friend class NewLlhTH2FFCN;

  NewLlhTH2F();
  virtual ~NewLlhTH2F();

  void PrepareAnalysis();

  // These are exact copies from NewLlhEnergy, but are required because fcn_ is of
  //  a different type of function
  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double nSrc, double gamma);

  virtual double EvaluateLlh(double* parArray) {
    return EvaluateLlh(parArray[0], parArray[1]);
  }

  // New for TH2F source hypothesis:
  /* @brief Loads pre-convolved source PDFs (i.e. PSF convolved with signal 
  / distribution) and the Detector Respose histos
  */
  void LoadTH2Fpdfs(char *filenamePDF);

  vector<double> GetSpaceRatios() { return spaceRatioVect_;}
 

};

class NewLlhTH2F_ParTranslator : public ParTranslator {
 protected:
  vector<double> setWeightVect_; // only one entry for each dataset now
/*
  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;
*/

 public:
  virtual ~NewLlhTH2F_ParTranslator() { }

  virtual const vector<double>
    Translate(int llhIndex, const vector<double>& parGlobal) const;

  virtual void SetTranslator(MultiAnalysisSet* mASet) {
    SetTranslator(mASet->GetAnalysisSetVect());
  }
  virtual void SetTranslator(const vector<AnalysisSet*>& aSetVect);
};


#endif // LLH_NEWLLHTH2F_H_
