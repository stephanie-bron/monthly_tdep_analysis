#include "llhTimeDep/public/NewLlhGausTime.h"
#include "llhTimeDep/public/TimePdfCollection.h"

#include "TGraph.h"

#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

  // Welcome to LlhGausTime.C!
  
  // This is meant to take a series of (already loaded)
  // I3Events, and a signal location and tests for 
  // compatability with a single Gaussian flare (in time).
  
  // If you want to test for a flare from a specific time
  // and duration you can set it by hand in the macro you
  // use to set up the search:
  
  //  llhEnergyFn.SetParDef(2, mean, 54561.,54971., 5e6., true, false);
  //  llhEnergyFn.SetParDef(3, log10(sigma), 0., -8., 4., true, false);
  //
  //         ---> SetParDef(par#, init, step_size min, max, optFix, optAuto)
  
  // will look for a flare specifically at mean, with width sigma
  // and will not look for the best flare (optFix==true)


double NewLlhGausTimeFCN::operator() (const vector<double>& par) const
{
  double srcFrac = par[0]/ptr->nEventsTot_;   // relative weight of src term
  //cout << "par[0] " << par[0] << " ptr->nEventsTot_ " << ptr->nEventsTot_ << endl;
  //cout << "par[1] " << par[1] << endl;
  //cout << "par[2] " << par[2] << endl;
  //cout << "par[3] " << par[3] << endl;
  double f=0.;
  
  if (srcFrac > ptr->srcFracMax_) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,ptr->srcFracMax_);
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (ptr->optStoreRatios_) { ptr->eventRatioVect_.clear(); }

  double gamma = par[1];  // Spectral index
  //const Event* eTemp = ptr->selectedList_.GetEvent(0);
 
  //if (ptr->optUseEnergy_ && (gamma >= 4. || gamma < 1. || TMath::IsNaN(gamma) ) ) { 
  if (ptr->optUseEnergy_ && (gamma >= ptr->eProbVect_[0]->GetGammaMax() || gamma < ptr->eProbVect_[0]->GetGammaMin() || TMath::IsNaN(gamma) ) ) {
    f = 1e50;
    return f;
  }
  
  double tmin = ptr->tmin_;
  double tmax = ptr->tmax_;

  double LogLambda = 0.;
  double TimeRatio = 1.;

  double mean = par[2];    // The peak of the Gaussian, in MJD
  double sigma = pow(10.,par[3]); // I use the log10 of the sigma, since it can vary over many
                           // orders of magnitude.
    
  if ( TMath::IsNaN(mean) || TMath::IsNaN(sigma) ) { // Sometimes something causes the FCN
                                                     // to test crazy values, here we just
                                                     // say to give a very high answer.
    //cout << "Sigma: " << sigma << " = exp(" << par[3] << ") is causing troubles?" << endl;
    f = 1e50;
    return f;
  }
  
  
  if ( ptr->monitorLevel_ > 4 ) { cout << "new GaussianTimePdf( " << tmin << ", " << tmax << ", " << mean << ", " << sigma << ",1.)" << endl; }
  TimePdf * timePdf = new GaussianTimePdf(tmin, tmax, mean, sigma,1.);
  timePdf->SetLivetime( tmax-tmin );
  //vector<double> trange = ptr->i3Set_->GetEventTimeModulePtr()->GetTimeRange();
  //timePdf->CheckTimeBounds(trange[0], trange[1]);

  ptr->nevs_=0; // number of events which contribute S/B > 1

  double lh;
  //double b=1.;
  
  for (int i=0;i<ptr->nEventsSelected_;i++)  { // loop over selected events

    const Event* event = ptr->selectedList_.GetEvent(i);

    TimeRatio = timePdf->GetPdfValue( event->GetTime().GetMJD() ) * timePdf->GetLivetime();
 
    // Perhaps we can speed things up... if the TimeRatio==0, skip the rest of the calculations.
    // include an out for high monitoring levels...
    if (ptr->monitorLevel_ < 4 && TimeRatio < 1e-50) {
      LogLambda += log(1.-srcFrac);
      continue;
    }


    // first comes the event's energy term
    double eRatio=1.;
    if (ptr->optUseEnergy_) {
      eRatio = ptr->eProbVect_[i]->GetEnergyProbGamma(*event, gamma) /
	           ptr->eProbVect_[i]->GetEnergyProbBkg(*event);
    }
        
    //As implemented, local coord background terms will automagically be included in spaceRatoVect   
    lh = ( srcFrac * ( ptr->spaceRatioVect_[i] * eRatio * TimeRatio - 1. ) + 1. );
    //cout << "srcFrac " << srcFrac << " TimeRatio " << TimeRatio << " lh " << lh << endl;
    
    if(lh > 1) {ptr->nevs_++;}
    LogLambda += log( lh );

    if ( (ptr->monitorLevel_ > 3 && lh > 1.) || ptr->monitorLevel_ > 4 ) { 
                                      // useful monitoring but with massive event-by-event spew
                                      // (~200 lines) going through each FCN call (~100-300)
      const I3Event* i3ev = (dynamic_cast<const I3Event*>(event));
                                      
      cout << lh << " " << LogLambda << " " << ptr->spaceRatioVect_[i] << " " << 
              eRatio << " " << event->GetTime().GetMJD() << " " << 
              timePdf->GetPdfValue(event->GetTime().GetMJD()) << " " << 
              i3ev->GetParams().runID << " " << i3ev->GetParams().eventID << endl;
    }

    if (ptr->optStoreRatios_) {
      ptr->eventRatioVect_.push_back(ptr->spaceRatioVect_[i]*eRatio*TimeRatio);
    }

  } //end loop over events
  
  if (ptr->optimizeTolerance_ > 0.) { // correction for events skipped in optimization
    LogLambda += (ptr->nEventsTot_-ptr->nEventsSelected_)*log(1.-srcFrac);
  }

  if ( ptr->JimsTerm_ && ( (2.0 * sigma) < timePdf->GetLivetime() ) ) { 
    LogLambda += log( sigma * 2.0/ timePdf->GetLivetime() );
  }
  ptr->margValue_ = log( sigma * 2.0/ timePdf->GetLivetime() );
  // Jim's Spectral index penalty
  // We could expect a gamma between 2.0 and 2.7
  // he applies a Gaussian with sigma 0.2 units in spectrum outside of that range
  if (ptr->SpectralPenalty_ && ptr->optUseEnergy_) { 
    if (gamma < 2.) LogLambda -= (2.-gamma)*(2.-gamma)/0.08;
    else if (gamma > 2.7) LogLambda -= (2.7-gamma)*(2.7-gamma)/0.08;
  }

  f = - LogLambda;    // What Minuit minimizes: -log L

  if ( ptr->monitorLevel_>1 ) { // more useful monitoring, with only
                           // one printout per FCN call. Still lots of spew
                           // but it's pretty managable for one or a few trials.
    printf("LogLambda=%12.6lg : nSrc=%9.4lg : gamma=%5.4lg : par1=%8.8lg : par2=%8.4lg : nEvs=%i : fact=%g\n",
//	   LogLambda,par[0],par[1],par[2],pow(10.,par[3]),ptr->nevs_,log( sigma * 2.0/timePdf->GetLivetime() ) );
	   LogLambda, srcFrac*ptr->nEventsTot_, gamma, mean, sigma, ptr->nevs_, ptr->margValue_ );
  }  
    
  if (timePdf) { delete timePdf; } // otherwise we'd get a memory leak
  return f;
}


NewLlhGausTime::NewLlhGausTime() :
  i3Set_(NULL),
  useLCBkgProb_(false),
  nullTestStat_(NULL),
  pvalHisto_(NULL),
  //fitfn_(NULL),
  histoForProb_(false),
  optUseEnergy_(true),
  eMaxRatioWarnStatus_(-1),  // default -1 means warn every time
  optStoreRatios_(false),
  icstatWarnLevel_(0),
  nSrcMin_(0.),
  srcFracMax_(0.5), // default is half of total nEvents
  gammaMin_(0.),
  gammaMax_(0.),
  logLambdaBest_(0.),
  nSrcBest_(0.),
  gammaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  nEventsTot_(0),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.),
  timePdfType_(0),
  fitfrac_(0.001),
  seedWtMin(1000)
{
  fcn_ = new NewLlhGausTimeFCN();
  fcn_->Point(this);

  minuit_ = new TFitterMinuit();
  minuit_->SetMinuitFCN(fcn_);
  // minuit_ takes over fcn_, so later we only delete minuit_, not fcn_

  // Call Migrad with 500 iterations maximum 
  minuit_->SetMaxIterations(500);
  // Set error Definition (1 for Chi square; 0.5 for negative log likelihood) 
  minuit_->SetErrorDef(0.5); 
  // Set Print Level (-1 no output; 1 standard output)
  minuit_->SetPrintLevel(-1);

  optParAuto_[0] = true;
  optParAuto_[1] = true;
  optParAuto_[2] = true;
  optParAuto_[3] = true;

  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}


void NewLlhGausTime::OptimizeEventSelection() {
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
    const I3Event* event = (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
    
    if (!event) { log_fatal("OptimizeEventSelection: failed cast to event.\n");}

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
        
    //cout << sigSpaceProb << "/" << bkgSpaceProb << " is it the LCBkgProb?" << endl;
    
    if ( useLCBkgProb_ ) {
        double lcProb;
        if (UseFolded_){
            lcProb = lcBkgProb_->BackgroundLCProb_folded(*event);
        }  else {
            lcProb = lcBkgProb_->BackgroundLCProb(*event); 
        }
        bkgSpaceProb = bkgSpaceProb*lcProb;
    }
    //cout << "nope..." << endl;    
    
    //    if event is in region with b=0 (what to do?)
    if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
    double spaceRatio = sigSpaceProb / bkgSpaceProb;

    double eMaxRatio;
    const EnergyProb* eProb(NULL);
    if (optUseEnergy_) { 
      eProb = event->GetEnergyProbFn();
      eMaxRatio = eProb->GetEnergyMaxRatio(*event);
      if (eMaxRatio <= 0.) {
	// Decide what sort of warning to give:
	if (eMaxRatioWarnStatus_ < 0) {          // Always warn
	  log_error("Error: gamma table is zero for this energy.  "
		    "Max Ratio=0.\n");
	  cout << event->GetParams().energyValue << endl;
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
      eVect_.push_back(*event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
      if (monitorLevel_ > 2) {
        cout << "(" << eventRaDeg << "," << eventDecDeg << ")  Sr : " << sigSpaceProb << "/" << bkgSpaceProb << " eP : " << eMaxRatio << endl;
      }
    }
  }

  if (monitorLevel_ > 0) {
    printf("Optimizing: %d events selected with maxRatio >%lg"
	   "  out of %d Ntotal.\n",
	   selectedList_.GetSize(), threshold, nEventsTot_);
  }
}



void NewLlhGausTime::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eVect_.clear();
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
      double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
      //    if event is in region with b=0 (what to do?)
      if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
      double spaceRatio = sigSpaceProb / bkgSpaceProb;
      const EnergyProb* eProb(NULL);
      if (optUseEnergy_) { 
	eProb = event->GetEnergyProbFn();
      }
      eVect_.push_back(*event);
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
  }
  
  nEventsSelected_ = eVect_.size();
  
}


void NewLlhGausTime::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  timePdfType_=1; //I should get rid of this, but set it here for now...

  PrepareAnalysis();
  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:
  assert(eVect_.size());
  //assert(selectedList_.GetSize());


  // sigmamin_ is set up to be the minimimum witdh of the
  // time Gaussian. If you let the width become arbitrarily small
  // the method could find that it can keep shrinking onto one strong event
  // and keep getting a better test statistic, which isn't what we want.
  // Sets a lower bound of the rms between the closest two events in time.

  sigmamin_ = 1000.;
  double timed;
  int sp = 1;
  vector<Time> ttt;

  for (int i=0;i<nEventsSelected_;i++) {
    ttt.push_back( eVect_[i].GetTime() );
  }

  for (unsigned int i=sp; i<ttt.size(); i++) {//This assumes events are time-ordered.
    timed = ttt[i].GetTime() - ttt[i-sp].GetTime();
    if (monitorLevel_ > 2) { cout << ttt[i].GetTime() << " "; }
    if (fabs(timed) < sigmamin_) {sigmamin_ = fabs(timed);}
  }

  sigmamin_ /= sqrt(2.);
  if (monitorLevel_ > 2) { cout << " sigmamin= " << sigmamin_ << endl; }
  ttt.clear();

  // Wipes out parameters (minuit crashes easily when changing existing pars)
  minuit_->Clear();
  minuit_->CreateMinimizer();  // default is kMigrad


  // DEFINE PARAMETERS

  GetFlareGuess(optUseEnergy_,1, nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_);

  if (optParAuto_[0]) {
    double nSrcMax = srcFracMax_*nEventsTot_;
    nSrcMin_ = 0;
    if (!nsrcGuess_) { nsrcGuess_=2.5; }
    minuit_->SetParameter(0, "nSrc", nsrcGuess_, 0.1, nSrcMin_, nSrcMax);
    //cout << "par set " << 0<< " " << "nSrc"<< " " << nsrcGuess_<< " " << 0.1<< " " <<nSrcMin_<< " " << nSrcMax << endl;
  }
 
  if (optParAuto_[1]) {
    if (optUseEnergy_) {
      gammaMin_ = eProbVect_[0]->GetGammaMin();
      gammaMax_ = eProbVect_[0]->GetGammaMax();
     
      //gammaInit = gammaGuess_;
      minuit_->SetParameter(1, "gamma", gammaGuess_, 0.1, gammaMin_, gammaMax_);
      //cout << "par set " << 1<< " " << "gamma"<< " " << gammaGuess_<< " " << 0.1<< " " <<gammaMin_<< " " << gammaMax_ << endl;

    } else {
      // If not using energy, then FCN will actually replace the
      // energy prob term internally... so these values don't really matter
      // Just fix the gamma so that Minuit doesn't try to minimize w.r.t. it
      minuit_->SetParameter(1, "gamma", 0., 0., 0., 0.);
    }
  }
  
  if (optParAuto_[2]) {
      double lowLimit = tmin_;
      double upLimit =  tmax_;
      double  initStepSize = sigmaGuess_;
      double  initValue = meanGuess_;
      minuit_->SetParameter(2, "mean", initValue, initStepSize, lowLimit, upLimit);
      //cout << "par set " << 2<< " " << "mean"<< " " <<initValue << " " << initStepSize<< " " << lowLimit<< " " << upLimit << endl;

  }
   
  if (optParAuto_[3]) {
      // Using the log10 of sigma as the fit parameter to make things more smooth.
      double initValueS = log10(sigmaGuess_);
      double lowLimitS = log10(sigmamin_);
      double upLimitS = log10(livetime_/2.);
      double initStepSizeS = 1.;
      minuit_->SetParameter(3, "sigma", initValueS, initStepSizeS, lowLimitS, upLimitS);
      //cout << "par set " << 3<< " " << "sigma"<< " " << initValueS<< " " << initStepSizeS<< " " << lowLimitS<< " " << upLimitS << endl;
  }

  // MINIMIZE

  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  minuit_->Minimize();
  stopwatch_minuitMigrad_.Stop();

  // GET RESULTS

  StoreLogLambdaBest();   // Set logLambdaBest_
  nSrcBest_  = GetPar(0);
  gammaBest_ = GetPar(1);
  meanBest_  = GetPar(2);
  sigmaBest_ = pow( 10., GetPar(3) );

  /* Handled in StoreLogLambdaBest(), I hope...

  // fix, because minimizer may sometimes give slightly negative value instead
  // of exact zero, which will cause probability calcultion to choke...
  if (logLambdaBest_ < 0.) {
    // we can always do better, since LogLambda(ns=0,gamma) = 0.
    logLambdaBest_ = 0.;
    nSrcBest_ = 0.;
    gammaBest_ = 0.;  // i.e., no fit makes sense, if nSrc=0.
  }
  */

  chiSq_ = 2.*logLambdaBest_;

  double p_temp;
  
  chisq_prob(chiSq_, ndof_, &p_temp, &chiSqProb_);
  
  estProb_ = chiSqProb_ / 2.;  // one-sided chi-sq prob

  stopwatch_MaximizeLlh_.Stop();
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  END MAXZIMIZELLH()  ///////////////////////
///////////////////////////////////////////////////////////////////////////////

// This is clumsy but seems like the only way to get this info
void NewLlhGausTime::StoreLogLambdaBest() 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx;
  minuit_->GetStats(amin, edm, errdef, nvpar, nparx);
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}

// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

void NewLlhGausTime::GetFlareGuess(bool useE, int useSpace, double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms) 
{

  vector<double> tVectortemp;
  vector<double> tVectorclose;

  double llhMax= -100.;
  double llhtemp, rms, avgtime;
  for (int j=0; j<nEventsSelected_; j++) {      
    if ( (eProbVect_[j]->GetEnergyMaxRatio(eVect_[j])*spaceRatioVect_[j]) > seedWtMin ) {
        tVectorclose.push_back(eVect_[j].GetMJD());
    }
  }
  
  sort(tVectorclose.begin(),tVectorclose.end());

  if (monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". Sigmamin is " << sigmamin_ << endl; }


//  const int ngamma=2;
//  double g[ngamma] = {2.0, 3.5};
  const int ngamma=1;
  double g[ngamma] = {2.0};

  int kmax = (int) close_;
  if (11 < kmax) { kmax = 11; }

  int n[] =  {1,2,3,4,5,10,15,20,25,30,35,40};
  int ns[] = {1,2,3,4,5, 5, 5, 5, 5, 5, 5, 5};

  int kmin=0;
  
  while ( tVectorclose.size() < n[kmax] ) { kmax--; }
  
    //cout << "kmax " << kmax << endl;
  for (int k=kmin; k<kmax; k++){ //k==1 for pairs, k==2 for triples, etc.
  
    for (unsigned int i=0; i<(tVectorclose.size()-n[k]); i++) {
      for (int j=0;j<=n[k];j++){
        tVectortemp.push_back(tVectorclose[i+j]);
      }
            
   //get mean of the times to seed this gaussian
      double sum=0, sumsq=0;
      llhtemp = 0.;
      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sum += tVectortemp[ii];
      }
      avgtime = sum/tVectortemp.size();

      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
      }

    //calculate the rms to use as the sigma
      rms = sqrt( sumsq/tVectortemp.size() );
      //cout << " GUESS k " << k << " i " << i << " avgtime " << avgtime << " rms " << rms << endl;

      if (useE) {
        for (int h=0; h<ngamma; h++) {
          llhtemp = EvaluateLlh( ns[k]+1., g[h], avgtime, rms );
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_mean = avgtime;
            Guess_rms = rms;
            Guess_nsrc = ns[k]+1.;
          }
        }
        tVectortemp.clear();
      } else { // I have no idea why you'd want to do the analysis without energy, but here you go...
        llhtemp = EvaluateLlh( ns[k]+1., 0., avgtime, rms );
        if (llhtemp > llhMax) {
          llhMax = llhtemp;
          Guess_mean = avgtime;
          Guess_rms = rms;
          Guess_nsrc = ns[k]+1.;
        }
        tVectortemp.clear();
      }
    }
  }
  
  tVectorclose.clear();
  
  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, Guess_rms );
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }
  
  if (useE || useSpace) { } //nip compiler complaints
      //cout << " GUESS " << Guess_nsrc<< " " << Guess_gamma<< " " << Guess_mean<< " " << Guess_rms << endl;

}


void NewLlhGausTime::SetNullTestStat(TH1D * inputhisto) {
  
  // This reads in a specific TH1D as the null test statistic
  // to use for p-values instead of using a chisquare distribution
  // it also fits an exponential to the upper tail (default set to top 0.01%).
  
  // It takes in the raw test statistic distribution (pdf)
  // and then makes the cdf to work with later.
  
  histoForProb_ = true;
  
  pvalHisto_ = new TH1D();
  nullTestStat_ = new TH1D();
  
  char sc[] = "Scale";
  //double firstfit=0;
  pvalHisto_ = DescendingCumulate(inputhisto,sc);
  nullTestStat_ = (TH1*)inputhisto->Clone("hnew");
  
  //int bins = inputhisto->GetNbinsX();
  //for (int i=1;i<bins;i++) {
  //  if (pvalHisto->GetBinContent(i) < fitfrac_){
  //    firstfit = pvalHisto->GetBinCenter(i);
  //    break;
  //  }
  //}
   
  //double max = 2*inputhisto->GetBinLowEdge(bins) - inputhisto->GetBinLowEdge(bins-1);
  //fitfn_ = new TF1("fitfn_","exp([0]+[1]*x)",firstfit,100);
   
  //pvalHisto->Fit(fitfn_,"QO","",firstfit,max);
  
}

double NewLlhGausTime::GetProbFromHisto(double teststat){ //, bool useFit) {

  // this looks at a null teststat distribution
  // which we loaded, and if it's in the exponential
  // tail, should go with a fit expo->Eval(teststat) (top 1%)
  
  int bin = pvalHisto->FindBin(teststat);
  double ptemp = pvalHisto->GetBinContent(bin);
  
//  if (ptemp < fitfrac_ && useFit) { //small enough that we should do the fit stuff
//    ptemp = fitfn_->Eval(teststat);
//  }
  return ptemp;
  
  // in case we don't feel like using a fit and it is higher than the distribution
  // just return 1/trials (conservative?)
  if( !pvalHisto->GetBinContent(bin) ) { return 1.0/nullTestStat_->GetEntries(); }
  
  return -1.; //if we get here something is wrong, 
              // and may as well hear about it trying to take the log
}


double NewLlhGausTime::EvalFCN(const vector<double>& parVect) const {
  return (*fcn_)(parVect);
  
/*  assert(nEventsTot_);
  int npar = 4;
  double gin[4];
  double par[4];
  double f;
  int iflag = 0;
  par[0] = parVect[0];
  par[1] = parVect[1];
  par[2] = parVect[2];
  par[3] = parVect[3];
  EvalMinuitFCN(npar, gin, f, par, iflag);
  return (-f);  // We want Llh, not -Llh */
  
}



double NewLlhGausTime::EvaluateLlh(double nSrc, double gamma, double mean, double sigma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  parVect.push_back(mean);
  parVect.push_back( log10(sigma) );
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}


/*
void NewLlhGausTime::EvalMinuitFCN (int &npar, double *gin, double &f, 
			       double *par, int iflag)
{
  double srcFrac = par[0]/nEvents_;   // relative weight of src term
  
  if (srcFrac > srcFracMax_) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,srcFracMax_);
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (storeRatios_) { eventRatio_.clear(); }

  double gamma = par[1];  // Spectral index
 
  if (useEnergy_ && (gamma > 4. || gamma < 1. || TMath::IsNaN(gamma) ) ) { 
    f = 1e50;
    return;
  }

  double LogLambda = 0.;
  double TimeRatio = 1.;
  double sigma = 0.;
  
  TimePdf * tPdf = NULL;

  double mean = par[2];    // The peak of the Gaussian, in MJD
  sigma = pow(10.,par[3]); // I use the log10 of the sigma, since it can vary over many
                           // orders of magnitude.
    
  if ( TMath::IsNaN(mean) || TMath::IsNaN(sigma) ) { // Sometimes something causes the FCN
                                                     // to test crazy values, here we just
                                                     // say to give a very high answer.
                                                     
    //cout << "Sigma: " << sigma << " = exp(" << par[3] << ") is causing troubles?" << endl;
    f = 1e50;
    return;
  }
    
  tPdf = new GaussianTimePdf(tmin_,tmax_,mean,sigma,1.);
  tPdf->SetLivetime(tmax_-tmin_);

  nevs_=0; // number of events which contribute S/B > 1

  double lh;
  
  for (int i=0;i<nEventsSelected_;i++)  { // loop over selected events
  
    // first comes the event's energy term
  
    double eRatio;
    if (useEnergy_) {
      eRatio = 
         eProb_->GetEnergyProbGamma(eVect_[i],gamma) /
         eProb_->GetEnergyProbBkg(eVect_[i]);
    } else {
      eRatio = 1.;
    }
    if ( monitorLevel_ > 3 ) { cout << eRatio << endl; }
    
    // local coordinate background stuff
    
    double b;
    if (useBkgPdfFunc_) {
      b = i3Set_->GetBackgroundTimeFunc()->BackgroundPdfFunc(eVect_[i],tPdf->GetLivetime());
    } else { 
      b = BkgPdfFunc(eVect_[i]->GetTime(),tPdf->GetLivetime());
    }
    
    if ( monitorLevel_ > 3 ) { cout << b << endl; }
    
    // Calculate the event's time term
    
    if (b) { TimeRatio = tPdf->GetPdfValue(eVect_[i]->GetTime()) / b; }
    else { TimeRatio = 1.; }
       
    lh = ( srcFrac * ( spaceRatio_[i] * eRatio * TimeRatio - 1. ) + 1. );
    
    if(lh > 1) {nevs_++;}
    LogLambda += log( lh );

    if ( monitorLevel_ > 3 && lh > 1. ) { // useful monitoring but with massive event-by-event spew
                                          // (~200 lines) going through each FCN call (~100-300)
      cout << lh << " " << LogLambda << " " << srcFrac << " " << spaceRatio_[i] << " " << eRatio << " " << TimeRatio << " " << b << " " << eVect_[i].GetParams().runID << " " << eVect_[i].GetParams().eventID << endl;
    }

    if (storeRatios_) {
      eventRatio_.push_back(spaceRatio_[i]*eRatio*TimeRatio);
    }

  } //end loop over events
  
  if (optimizeTolerance_ > 0.) { // correction for events skipped in optimization
    LogLambda += (nEvents_-nEventsSelected_)*log(1.-srcFrac);
  }

  if ( timePdfType_ == 1 ) 
  { 
      if ( JimsTerm_ && ( 2.*sigma<tPdf->GetLivetime() ) ) { 
        LogLambda += log( 2.*sigma/tPdf->GetLivetime() ); 
      }
      margValue_ = log( sigma * 2.0/tPdf->GetLivetime() );
//      blocksThreshold_ = log( sigma * sqrt( 2.0*TMath::Pi())/tPdf->GetLivetime() );//par[0] ) );
  } 

  // Jim's Spectral index penalty
  // We could expect a gamma between 2.0 and 2.7
  // he applies a Gaussian with sigma 0.2 units in spectrum outside of that range
  if (SpectralPenalty_ && useEnergy_) { 
    if (gamma < 2.) LogLambda -= (2.-gamma)*(2.-gamma)/0.08;
    else if (gamma > 2.7) LogLambda -= (2.7-gamma)*(2.7-gamma)/0.08;
  }

  f = - LogLambda;    // What Minuit minimizes: -log L

  if ( monitorLevel_>1 ) { // more useful monitoring, with only
                           // one printout per FCN call. Still lots of spew
                           // but it's pretty managable for one or a few trials.  
    printf("LogLambda=%12.6lg : nSrc=%9.4lg : gamma=%5.4lg : par1=%8.5lg : par2=%8.4lg : nEvs=%i : fact=%g\n",
	   LogLambda,par[0],par[1],par[2],pow(10.,par[3]),nevs_,margValue_);
  }  
    

  if (npar || gin || iflag) { }; // touch these variables just to
  // eliminate warnings during compile
  
  if (tPdf) { delete tPdf; }

} //*/


// TGraph* NewLlhGausTime::GetContour(double sigma, int npoints, int pa1, int pa2) {
//   minuit_->SetFCN(MinuitWrapperFCN);
//   SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

//   minuit_->SetErrorDef(sigma*sigma);
//   TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

//   SetMinuitWrapperPtr(NULL);  // reset pointer
//   return g;
// }






// For MultiAnalysis: Translate Global Par values to Individual Par Values:
const vector<double> NewLlhGausTime_ParTranslator::Translate(int llhIndex, const vector<double>& parGlobal) const {
    double nSrcGlobal = parGlobal[0];
    double gamma = parGlobal[1];
    vector<double> parLocal(4);
    double weight =  LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);
    // recall: srcWeightVect_[i] is itself a vector<double> 
    parLocal[0] = nSrcGlobal * weight;
    parLocal[1] = gamma;
    parLocal[2] = parGlobal[2];
    parLocal[3] = parGlobal[3];
    return parLocal;
}
  

void NewLlhGausTime_ParTranslator::SetTranslator(const vector<AnalysisSet*>& aSetVect) {
  int nSets = aSetVect.size();

  srcWeightVect_.clear();
  vector<double> tempVect(nStopsGamma_,0);
  for (int i=0; i<nSets; ++i) {
    srcWeightVect_.push_back( tempVect );
  }

  for (int g = 0; g < nStopsGamma_; ++g) {
    double gamma = gammaMin_ + 
      (gammaMax_-gammaMin_) * (double(g))/(nStopsGamma_-1);
    double sum = 0.;
    for (int i=0; i<nSets; ++i) {
      double nev = 
	aSetVect[i]->GetMeanSrcNevForFluxModel(PowerLawFlux(1,-gamma));
      srcWeightVect_[i][g] = nev;
      sum += nev;
    }
    // for each choice of gamma, normalize weights of different data sets to 1
    cout << "Weights for Gamma=" << gamma << ": ";
    for (int i=0; i<nSets; ++i) {
      srcWeightVect_[i][g] /= sum;
      cout << srcWeightVect_[i][g] << " ";
    }
    cout << endl;
  }    
}


