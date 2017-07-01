#include "llh/public/LlhEnergy.h"

#include "TMinuit.h"
#include "TGraph.h"

#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

LlhEnergy::LlhEnergy() :
  AnalysisLlh(2),  // initialize minuit_ with 2 parameters for minimizing
  useEnergy_(true),
  eMaxRatioWarnStatus_(-1),  // default -1 means warn every time
  storeRatios_(false),
  icstatWarnLevel_(0),
  nSrcMin_(0.),
  srcFracMax_(0.5), // default is half of total nEvents
  gammaMin_(0.),
  gammaMax_(0.),
  logLambdaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  nSrcBest_(0.),
  gammaBest_(0.),
  nEventsTot_(0),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.)
{
  parDefArray_[0].name = "nSrc";
  parDefArray_[1].name = "gamma";

  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}


void LlhEnergy::OptimizeEventSelection() {
  const EquatorialDeg *srcEqDeg = 
    ( dynamic_cast<const EquatorialDeg*>(srcCoord_) );

  double decMinDeg = -90.;
  double decMaxDeg =  90.;
  double raRangeDeg = 360.;
  if (optimizeAngleDeg_ > 0.) {
    decMinDeg = srcEqDeg->GetDec() - optimizeAngleDeg_;
    decMaxDeg = srcEqDeg->GetDec() + optimizeAngleDeg_;
    if (decMinDeg < -90.) { decMinDeg = -90.; }
    if (decMaxDeg >  90.) { decMaxDeg =  90.; }
    // scale the ra range according to the dec closest to one of the poles
    double cosFactor = cos(decMinDeg*TMath::DegToRad());
    double cosFactorAlt = cos(decMaxDeg*TMath::DegToRad());
    if (cosFactorAlt < cosFactor) { cosFactor = cosFactorAlt; }
    if (cosFactor>0) {
      raRangeDeg = optimizeAngleDeg_ / cosFactor;
    }
  }

  double threshold = 
    (1./srcFracMax_ - 1.) * ( exp(optimizeTolerance_/nEventsTot_) -1. );

  const EventPtrList* evList = aSet_->GetEventPtrList();

  for (int i=0; i<nEventsTot_; ++i) {
    const I3Event* event = 
      (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
    assert(event);

    double eventRaDeg = event->GetEquatorialDeg().GetRa();
    // use this function to test whether ra is within range
    if ( ModDistanceFn(eventRaDeg,srcEqDeg->GetRa(),360.) > raRangeDeg ) {
      continue;  // skip this event
    }

    double eventDecDeg = event->GetEquatorialDeg().GetDec();
    if (eventDecDeg < decMinDeg || eventDecDeg > decMaxDeg) {
      continue;  // skip this event
    }

    double sigSpaceProb = event->ProbFrom(*srcCoord_);
    double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
    assert(bkgSpaceProb); // if event is in region with b=0 (what to do?)
    double spaceRatio = sigSpaceProb / bkgSpaceProb;

    double eMaxRatio;
    const EnergyProb* eProb(NULL);
    if (useEnergy_) { 
      eProb = event->GetEnergyProbFn();
      eMaxRatio = eProb->GetEnergyMaxRatio(*event);
      if (eMaxRatio <= 0.) {
	// Decide what sort of warning to give:
	if (eMaxRatioWarnStatus_ < 0) {          // Always warn
	  log_error("Error: gamma table is zero for this energy.  "
		    "Max Ratio=0.\n");
	} else if (eMaxRatioWarnStatus_ == 0) {  // Warn Once
	  log_error("***\nWARNING!\nAt LEAST one data event has zero for "
		    "its gamma energy pdf.\nAll such events will be "
		    "effectively ELIMINATED when using llh with energy.\n"
		    "Llh is currently set so that you will receive this "
		    "warning ONLY ONE TIME.\n***\n");
	  eMaxRatioWarnStatus_ = 1;
	}
	// else don't warn again
      }
    } else {
      eMaxRatio = 1.;
    }
    
    if (spaceRatio*eMaxRatio > threshold) {  // select this event
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
    
    if (monitorLevel_ > 2) {
      cout << "(" << eventRaDeg << "," << eventDecDeg << ")";
      cout << "  Sr : " << sigSpaceProb << "/" << bkgSpaceProb;
      cout << " eP : " << eProb << endl;
    }
  }

  if (monitorLevel_ > 0) {
    printf("Optimizing: %d events selected with maxRatio >%lg"
	   "  out of %d Ntotal.\n",
	   selectedList_.GetSize(), threshold, nEventsTot_);
  }
}



void LlhEnergy::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  // LOAD DATA

  assert(aSet_);
  assert(srcCoord_);

  const EventPtrList* evList = aSet_->GetEventPtrList();
  assert(evList);
  nEventsTot_ = evList->GetSize();
  assert(nEventsTot_);  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eProbVect_.clear();
  spaceRatioVect_.clear();

  eventRatioVect_.clear();


  if (optimizeTolerance_ > 0.) {
    stopwatch_optimize_.Start(false);  //  false = don't reset
    OptimizeEventSelection();
    stopwatch_optimize_.Stop();
  } 
  else 
  { // no optimization
    for (int i=0; i<nEventsTot_; ++i) {
      const I3Event* event = 
	(dynamic_cast<const I3Event*> (evList->GetEvent(i)));
      assert(event);
      double sigSpaceProb = event->ProbFrom(*srcCoord_);
      double bkgSpaceProb = 
	event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
      assert(bkgSpaceProb); // if event is in region with b=0 (what to do?)
      double spaceRatio = sigSpaceProb / bkgSpaceProb;
      const EnergyProb* eProb(NULL);
      if (useEnergy_) { 
	eProb = event->GetEnergyProbFn();
      }
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
  }

  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:

  assert(selectedList_.GetSize());

  // DEFINE PARAMETERS

  if (parDefArray_[0].optAuto ) {
    double nSrcMax = srcFracMax_*nEventsTot_;
    nSrcMin_ = 0.;
    // SetParDef(i, initVal, initStepSize, lowLim, upLim, optFix, optAuto);
    SetParDef(0, nSrcMax/2., 0.1, nSrcMin_, nSrcMax, false, true);
  }

  if (parDefArray_[1].optAuto ) {
    if (useEnergy_) {
      // Use eProb of *first* event... assume gamma range same for all????
      gammaMin_ = eProbVect_[0]->GetGammaMin();
      gammaMax_ = eProbVect_[0]->GetGammaMax();
      double gammaInit = (gammaMin_+gammaMax_)/2.;

      // SetParDef(i, initVal, initStepSize, lowLim, upLim, optFix, optAuto);
      SetParDef(1, gammaInit, 0.1, gammaMin_, gammaMax_, false, true);
    } else {
      // If not using energy, then FCN will actually replace the
      // energy prob term internally... so these values don't really matter
      // Just fix the gamma so that Minuit doesn't try to minimize w.r.t. it
      SetParDef(1, 0., 0., 0., 0., true, true);
    }
  }

  AssignMinuitParDef();

  // MINIMIZE

  // Set the Minuit FCN indirectly, by using a wrapper fcn which 
  // internally will point to 'this' object and its own EvalMinuitFCN

  minuit_->SetFCN(MinuitWrapperFCN);
  SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN
  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  minuit_->Migrad();          // do the minimization
  stopwatch_minuitMigrad_.Stop();
  SetMinuitWrapperPtr(NULL);  // reset pointer


  // GET RESULTS

  double fmin, fedm, errdef;
  int nvpar, nparx, icstat;
  minuit_->mnstat(fmin, fedm, errdef, nvpar, nparx, icstat);

  if (icstat <= icstatWarnLevel_) {
    log_warn("icstat=%d : ",icstat);
    if (icstat == 0) { log_warn("covar. matrix not calculated at all.\n"); }
    if (icstat == 1) { log_warn("covar. matrix is approx. only\n"); }
    if (icstat == 2) { log_warn("full covar.matrix, but forced pos-def.\n");}
    if (icstat == 3) { log_warn("full, accurate covariance matrix.\n");}
  }

  logLambdaBest_ = -fmin;
  nSrcBest_ = GetPar(0);
  gammaBest_ = GetPar(1);

  if (logLambdaBest_ < 0.) {
    // we can always do better, since LogLambda(ns=0,gamma) = 0.
    nSrcBest_ = 0.;
    gammaBest_ = 0.;  // i.e., no fit makes sense, if nSrc=0.
    logLambdaBest_ = 0.;
  }

  chiSq_ = 2.*logLambdaBest_;

  double p_temp;
  if (useEnergy_) {
    chisq_prob(chiSq_, 2., &p_temp, &chiSqProb_);
  } else {
    chisq_prob(chiSq_, 1., &p_temp, &chiSqProb_);
  }


  stopwatch_MaximizeLlh_.Stop();
}


TGraph* LlhEnergy::GetContour(double sigma, int npoints, int pa1, int pa2) {
  minuit_->SetFCN(MinuitWrapperFCN);
  SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

  minuit_->SetErrorDef(sigma*sigma);
  TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

  SetMinuitWrapperPtr(NULL);  // reset pointer
  return g;
}


double LlhEnergy::EvaluateLlh(double ns, double gamma) {
  assert(nEventsTot_);
  int npar = 2;
  double gin[2];
  double f;
  double par[2];
  int iflag = 0;
  par[0] = ns;
  par[1] = gamma;
  storeRatios_ = true;
  EvalMinuitFCN(npar, gin, f, par, iflag);
  storeRatios_ = false;
  return (-f);  // We want Llh, not -Llh
}



void LlhEnergy::EvalMinuitFCN (int &npar, double *gin, double &f, 
			       double *par, int iflag)
{
  double srcFrac = par[0]/nEventsTot_;   // relative weight of src term
  if (srcFrac > srcFracMax_) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,srcFracMax_);
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (storeRatios_) { eventRatioVect_.clear(); }

  double gamma = par[1];                // center of Gaussian

  double LogLambda=0.;

  for (int i=0; i<selectedList_.GetSize(); ++i) {
    
    double eRatio;
    if (useEnergy_) {
      const Event* event = selectedList_.GetEvent(i);
      eRatio = eProbVect_[i]->GetEnergyProbGamma(*event, gamma) /
	eProbVect_[i]->GetEnergyProbBkg(*event);
    } else {
      eRatio = 1.;
    }
    // CAN USE THIS STOPWATCH TO ISOLATE ONE PIECE OF THE LLH CALC.
    //  stopwatch_minuitFCN_.Start(false);  //  false = don't reset
  
    LogLambda += log( srcFrac * ( spaceRatioVect_[i] * eRatio - 1) + 1);

    // DON'T FORGET TO STOP AFTERWARDS!
    //    stopwatch_minuitFCN_.Stop();

    if (storeRatios_) {
      eventRatioVect_.push_back(spaceRatioVect_[i]*eRatio);
    }      
  }

  if (optimizeTolerance_ > 0.) { // correction for events skipped in optimization
    LogLambda += (nEventsTot_ - selectedList_.GetSize())*log(1.-srcFrac);
  }

  f = - LogLambda;    // What Minuit minimizes: -log L

  if (monitorLevel_>1) {
    printf("LogLambda=%12.6lg  :  nSrc=%9.4lg  :  gamma=%5.4lg\n",
	   LogLambda,par[0],par[1]);
  }

  if (npar || gin || iflag) { }; // touch these variables just to
  // eliminate warnings during compile

}
