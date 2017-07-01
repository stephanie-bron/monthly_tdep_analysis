#include "llh/public/NewLlhTH2F.h"

#include "TGraph.h"
#include "TFile.h"

#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

#include "llh/public/CoordinateTransform.h"


NewLlhTH2F::NewLlhTH2F()
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
  optimizeAngleDeg_ = 0.;

  //fcn_ = new NewLlhTH2FFCN();
  fcn_ = new NewLlhEnergyFCN();
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

NewLlhTH2F::
~NewLlhTH2F()
{
}

void NewLlhTH2F::OptimizeEventSelection() {
  log_fatal("Optimiziation not implemented for stacking yet!\n");
  return;
}


void NewLlhTH2F::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eProbVect_.clear();
  spaceRatioVect_.clear();
  eventRatioVect_.clear(); //this gets filled in NewLlhEnergy.C


  optimizeTolerance_ = 0;
  // no optimization
  double evra=0, evdec=0;
  double evX=0, evY=0;
  double sigmaMin = 0; // Hardcoded values for now...
  double sigmaMax = 3;
  const int sigmaNbins_ = 30;

  for (int i=0; i<nEventsTot_; ++i) {
    const I3Event* event = (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
    assert(event);

    EquatorialDeg ed = event->GetEquatorialDeg();

    evra = ed.GetRa();
    evdec= ed.GetDec();
    evX = evra;  //everything in equatorial in degrees now
    evY = evdec; //everything in equatorial in degrees now

    double sigmaEv = (double)event->GetParams().parafitSigmaDeg;

    // Find which conv TH2F is appropriate for this sigma
    // Hardcoded binning for now
    int sigmaBin = (int)lround( (sigmaEv - sigmaMin - (sigmaMax-sigmaMin)/(2*sigmaNbins_)) * sigmaNbins_/(sigmaMax-sigmaMin) );
    if (sigmaBin > 29) sigmaBin = 29;

    //find right bin in the found histo
    int xIndex = hPdf_[sigmaBin].GetXaxis()->FindBin(evX);
    int yIndex = hPdf_[sigmaBin].GetYaxis()->FindBin(evY);

    //now space prob for event
    double sigSpaceProb = hPdf_[sigmaBin].GetBinContent(xIndex, yIndex); 
    double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
    //    if event is in region with b=0 (what to do?)
    if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
    
    //finally we have the spacial prob ratio for this event
    double spaceRatio = sigSpaceProb / bkgSpaceProb;

    const EnergyProb* eProb(NULL);
    if (optUseEnergy_) {
  	eProb = event->GetEnergyProbFn();
    }

    //insert event in list, space and energy probability in list
    selectedList_.AddEvent(event);
    eProbVect_.push_back(eProb);
    spaceRatioVect_.push_back(spaceRatio);
  } // end loop over events
}


double NewLlhTH2F::EvalFCN(const vector<double>& parVect) const {
  return (*fcn_)(parVect);
}

double NewLlhTH2F::EvaluateLlh(double nSrc, double gamma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}

void NewLlhTH2F::LoadTH2Fpdfs(char *filenamePDF) {

  // File containing pre-convolved TH2Fs
  TFile *f1 = new TFile(filenamePDF,"READ");
  const int sigmaNbins_ = 30;
  char nameBuf[30];
  for (int i=0; i<sigmaNbins_; i++) {
    sprintf(nameBuf,"h%d",i);
    string sBuf = nameBuf;
    // dynamic_cast the TObject (from TFile::Get) as pointer to TH2F, then de-ref
    hPdf_[i] = *(dynamic_cast<TH2F*>( f1->Get(sBuf.c_str()) ));
    hPdf_[i].SetDirectory(0); // The histogram needs to get decoupled from the ROOT file!
    //hPdf_[i].Scale(180.*180.*4/3.14);
    //hPdf_[i].Scale(1./hPdf_[i].Integral());
    //hPdf_[i].Scale(3.141592/hPdf_[i].Integral());
    cout<<sBuf<<" integrates to "<<hPdf_[i].Integral()<<endl;
  }
  f1->Close();
}

const vector<double> NewLlhTH2F_ParTranslator::Translate(int llhIndex, const vector<double>& parGlobal) const {
  double nSrcGlobal = parGlobal[0];
  double gamma = parGlobal[1];

  vector<double> parLocal(2);

  double weight = setWeightVect_[llhIndex];
  // recall: srcWeightVect_[i] is itself a vector<double> 

  // nSrc , must be scaled to relative weight of this source
  parLocal[0] = nSrcGlobal * weight;
  parLocal[1] = gamma;
  return parLocal;
}

void NewLlhTH2F_ParTranslator::SetTranslator(const vector<AnalysisSet*>& aSetVect) {
  int nSets = aSetVect.size();

  setWeightVect_.clear();
  setWeightVect_.resize(nSets);

  if (nSets>1){
    double sum = 0.;
    for (int i=0; i<nSets; ++i) {
      double nev = aSetVect[i]->GetMeanSrcNev(); // Assuming a fixed flux model!
      setWeightVect_[i] = nev;
      sum += nev;
    }
    // for each choice of gamma, normalize weights of different data sets to 1
    cout << "Weights for Datasets = " << flush;
    for (int i=0; i<nSets; ++i) {
      setWeightVect_[i] /= sum;
      cout << setWeightVect_[i] << " " << flush;
    }
  } else {
      setWeightVect_[0] = 1.0;
  }
  cout << endl;
}


