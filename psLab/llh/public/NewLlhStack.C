#include "llh/public/NewLlhStack.h"

#include "TGraph.h"

#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"




double NewLlhStackFCN::operator() (const vector<double>& par) const {
  double nSrc = par[0];
  double gamma = par[1];

  if (ptr->optStoreRatios_) { ptr->eventRatioVect_.clear(); }

  double srcFrac = nSrc/ptr->nEventsTot_;   // relative weight of src term
  if (srcFrac > ptr->srcFracMax_) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,ptr->srcFracMax_);
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  ///
  // JD: I use the convention that flux is proportional to E^gamma
  //  SimpleEnergyProb uses E^-gamma... (sorry for confusion)
  ///

  // round a double to nearest 'long'  This avoids round returning a double of
  //  0.99... which (int) might cast down to zero.
  int gammaIndex=(int)lround(10*(-1.*gamma+4)-0.5); // for look-up in Stacking table
  // For now, table size and binning is fixed in macro_StackedWeights.C
  if (gammaIndex<0) gammaIndex=0;
  if (gammaIndex>29) gammaIndex=29;

  double stackGammaMin = -4;
  double stackGammaMax = -1;
  int nGammaBins = 30; // find a global way of setting this...
  // This gamma index corresponds to a bin center of:
  double gammaCenterThisBin = (gammaIndex+0.5)*(stackGammaMax-stackGammaMin)/nGammaBins + stackGammaMin;


  double logLambda=0.;
  double ratioSum=0.;
  double srcWeightSum=0.;
  int srcIndex=0;

  for (int i=0; i<ptr->selectedList_.GetSize(); ++i) {
    double eRatio = 1.;
    if (ptr->optUseEnergy_) {
      const Event* event = ptr->selectedList_.GetEvent(i);
      eRatio = ptr->eProbVect_[i]->GetEnergyProbGamma(*event, gamma) /
	ptr->eProbVect_[i]->GetEnergyProbBkg(*event);
    }
    // CAN USE THIS STOPWATCH TO ISOLATE ONE PIECE OF THE LLH CALC.
    //  stopwatch_minuitFCN_.Start(false);  //  false = don't reset
  
    // Iterate over each source, summing up signal terms and src weights
    ratioSum=0.;
    srcWeightSum=0.;
    srcIndex=0;
    for (unsigned j=0; j<ptr->spaceRatioVects_.size(); j++) {
      double srcWeightInterpValue = 0;
      bool interpDone = 0;      // Nominal value of srcWeightTable before interpolation
      double valueAtThisBin = ptr->srcWeightsTable_[srcIndex][gammaIndex];
      if ( gammaIndex==0 && -1.*gamma<gammaCenterThisBin ) {
        double lowEdge = -4.;
        srcWeightInterpValue = valueAtThisBin * (-1.*gamma-lowEdge)/(gammaCenterThisBin-lowEdge);
        interpDone = 1;
      }
      if ( gammaIndex == 29 && -1.*gamma>gammaCenterThisBin ){
        double hiEdge = -1;
        srcWeightInterpValue = valueAtThisBin * (-1.*gamma-hiEdge)/(hiEdge-gammaCenterThisBin);
        interpDone = 1;
      }

      if (!interpDone) {
        int dir = 1;
        if (-1*gamma<gammaCenterThisBin) { dir = -1; }
        double valueAtNextBin = ptr->srcWeightsTable_[srcIndex][gammaIndex+dir];
        double gammaCenterNextBin = (gammaIndex + dir + 0.5)*(stackGammaMax-stackGammaMin)/nGammaBins + stackGammaMin;

        double fractionOffset = fabs((-1*gamma-gammaCenterThisBin)/(gammaCenterNextBin-gammaCenterThisBin));
        srcWeightInterpValue = valueAtThisBin + (valueAtNextBin-valueAtThisBin)*fractionOffset;
      }
      // Interpolation turned out to be important for MINUIT (smooth llh landscape)
//cerr << "i: " << i << " 1: " << srcWeightInterpValue << " 2: " << ptr->spaceRatioVects_[srcIndex][i] << " 3: " << eRatio << endl;
assert(ptr->spaceRatioVects_[srcIndex][i]==ptr->spaceRatioVects_[srcIndex][i]);
      ratioSum += srcWeightInterpValue *
                    ptr->spaceRatioVects_[srcIndex][i]*eRatio;
      srcWeightSum += srcWeightInterpValue; // for normalization
      srcIndex++;
    }

//cerr << "1: " << srcFrac << " 2: " << ratioSum << " 3: " << srcWeightSum << endl;
    //logLambda += log( srcFrac * ( ptr->spaceRatioVect_[i] * eRatio - 1) + 1);
    logLambda += log( srcFrac * ( ratioSum/srcWeightSum - 1) + 1);

    // DON'T FORGET TO STOP AFTERWARDS!
    //    stopwatch_minuitFCN_.Stop();

    if (ptr->optStoreRatios_) {
      ptr->eventRatioVect_.push_back(ratioSum/srcWeightSum);
    }      
  }

  if (ptr->optimizeTolerance_ > 0.) { 
    // correction for events skipped in optimization
    log_fatal("This is skipped unless you implemented optimization!\n (in that case, double-check that this code still works)\n");
    logLambda += (ptr->nEventsTot_ - ptr->selectedList_.GetSize())*log(1.-srcFrac);
  }

  if (ptr->monitorLevel_>1) {
    printf("LogLambda=%12.6lg  :  nSrc=%9.4lg  :  gamma=%5.4lg\n",
	   logLambda, nSrc, gamma);
  }

  return -logLambda;   // What Minuit minimizes: -log L
}




NewLlhStack::NewLlhStack()
{
  optUseEnergy_ = true;
  eMaxRatioWarnStatus_ = -1;  // default -1 means warn every time
  optStoreRatios_ = false;
  icstatWarnLevel_ = 0;
  nSrcMin_ = 0.;
  srcFracMax_ = 0.5; // default is half of total nEvents
  gammaMin_ = 0.;
  gammaMax_ = 0.;
  logLambdaBest_ = 0.;
  nSrcBest_ = 0.;
  gammaBest_ = 0.;
  chiSq_ = 0.;
  chiSqProb_ = 0.;
  nEventsTot_ = 0;
  monitorLevel_ = 0;
  optimizeTolerance_ = 0.;
  optimizeAngleDeg_ = 0.;

  fcn_ = new NewLlhStackFCN();
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

  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}

void NewLlhStack::SetSourceCoords(vector<EquatorialDeg>& sourceCoords) {
  for ( vector<EquatorialDeg>::iterator iter = sourceCoords.begin();
    iter != sourceCoords.end(); iter++ )
  {
    srcCoords_.push_back(*iter);
  }
  return;
}

void NewLlhStack::SetSourceSigmas(vector<double>& sourceSigmas) {
  // In the past, this std::copy has caused seg faults...
  // Trying agin with resize first:
  srcSigmas_.resize(sourceSigmas.size());
  std::copy(sourceSigmas.begin(), sourceSigmas.end(), srcSigmas_.begin());
  return;
}

void NewLlhStack::SetStackedWeightTable(vector<vector<double> > srcWeightsArray) {
  srcWeightsTable_.clear();

 // This compiled but segfaulted:
 //copy(srcWeightsArray.begin(), srcWeightsArray.end(), srcWeightsTable_.begin());

  for ( vector<vector<double> >::iterator iter = srcWeightsArray.begin();
    iter != srcWeightsArray.end(); iter++ )
  {
    srcWeightsTable_.push_back(*iter);
  }

/*
  // Logging, if desired...
  for (vector< vector<double> >::size_type u = 0; u < srcWeightsTable_.size(); u++)
  {
      for (vector<double>::size_type v = 0; v < srcWeightsTable_[u].size(); v++) {
         cout << srcWeightsTable_[u][v] << " ";
      }
      cout << endl;
  }
*/
}



void NewLlhStack::OptimizeEventSelection() {
  log_fatal("Optimiziation not implemented for stacking yet!\n");
  return;
}


void NewLlhStack::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (srcCoords_.size() == 0) { log_fatal("PrepareAnalysis: srcCoords were not set.\n"); }
  if (srcSigmas_.size() == 0) { log_fatal("PrepareAnalysis: srcSigmas were not set.\n"); }
  if (srcWeightsTable_.size() == 0) { log_fatal("PrepareAnalysis: srcWeightsTable not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eProbVect_.clear();
  spaceRatioVect_.clear();
  spaceRatioVects_.resize( srcCoords_.size() );

  eventRatioVect_.clear();


  if (optimizeTolerance_ > 0.) {
    stopwatch_optimize_.Start(false);  //  false = don't reset
    log_error("Optimization not implemented yet! Setting Tolerance to 0.\n");
    optimizeTolerance_ = 0;
    stopwatch_optimize_.Stop();
  } 
  if (optimizeTolerance_ > 0.) {
    OptimizeEventSelection();
  }
  else
  { // no optimization
    int srcIndex=0;
    // Loop over sources
    for ( vector<EquatorialDeg>::iterator iter = srcCoords_.begin();
        iter != srcCoords_.end(); iter++ ) {
      spaceRatioVects_[srcIndex].clear();
      for (int i=0; i<nEventsTot_; ++i) {
        const I3Event* event =
    (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
        assert(event);

        // JD: I'm taking Spatial Prob out of "Event" and into "Llh"
        //  since it is a function of the size of source and the reco of event
        //double sigSpaceProb = event->ProbFrom(*iter); // OLD
        double r = event->GetCoord().DistanceTo(*iter);
        double es = event->GetParams().parafitSigmaDeg; // event sigma
        double ss = srcSigmas_[srcIndex]; // source sigma
        double sigma = sqrt( (es*es + ss*ss) ); // Convolved Gaussians
        double sigSpaceProb = exp(-r*r/(sigma*sigma*2)) / (2.*TMath::Pi()*sigma*sigma
);

        double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
        //    if event is in region with b=0 (what to do?)
        if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
        double spaceRatio = sigSpaceProb / bkgSpaceProb;
        spaceRatioVects_[srcIndex].push_back(spaceRatio);
      } // end loop over events
    srcIndex++;
    } // end loop over sources

    // Separate loop over events for energy and event list
    for (int i=0; i<nEventsTot_; ++i) {
      const I3Event* event =
  (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
      assert(event);
      const EnergyProb* eProb(NULL);
      if (optUseEnergy_) {
    eProb = event->GetEnergyProbFn();
      }
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
    }
  } // End of no optimization option

}

double NewLlhStack::EvalFCN(const vector<double>& parVect) const {
  return (*fcn_)(parVect);
}

double NewLlhStack::EvaluateLlh(double nSrc, double gamma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}
