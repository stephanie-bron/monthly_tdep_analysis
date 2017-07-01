
#include "llhTimeDep/public/MultiBlockAnalysisFn.h"
#include "iostream"

double MultiBlockAnalysisFCN::operator() (const vector<double>& parVect) const {
  double result = 0.;
  for (int i=0; i<int(ptr->analysisFnVect_.size()); ++i) {
    vector<double> individualParVect = ptr->parTrans_->Translate(i, parVect);
    result += ptr->analysisFnVect_[i]->EvalFCN(individualParVect);
  }
  return result;
}



MultiBlockAnalysisFn::MultiBlockAnalysisFn() {
  fcn_ = new MultiBlockAnalysisFCN();
  fcn_->Point(this);

  nPar = 4;
  laglimit = 0.5;
  Ndof=-1;
  minuit_ = new TFitterMinuit();
  minuit_->SetMinuitFCN(fcn_);
  // minuit_ takes over fcn_, so later we only delete minuit_, not fcn_

  // Call Migrad with 500 iterations maximum
  minuit_->SetMaxIterations(500);
  // Set error Definition (1 for Chi square; 0.5 for negative log likelihood)
  minuit_->SetErrorDef(0.5);
  // Set Print Level (-1 no output; 1 standard output)
  minuit_->SetPrintLevel(-1);
  // default is kMigrad
  minuit_->CreateMinimizer();
}


double MultiBlockAnalysisFn::EvalFCN(const vector<double>& parVect) const {
  return (*fcn_)(parVect);
}

double MultiBlockAnalysisFn::EvaluateLlh(double *parValueArray) {
  vector<double> parVect;

  //for (int i=0; i<int(parDefVect_.size()); ++i) {
  for (int i=0; i<nPar; ++i) {
    parVect.push_back(parValueArray[i]);
  }
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}


void MultiBlockAnalysisFn::StoreLogLambdaBest()
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



//    * Oh lordy, if we assume that we want to use something
//    * like the untriggered search over multiple years,
//    * that's a search which uses the event info to seed the
//    * first guess. I'm sure this will also be useful for plotting.
//    *
//    * This also assumes you don't want some sort of extra relative
//    * weighting that comes from the detector config/cuts, could be
//    * cool to have that too (in the events?).
//    *
/*
vector<I3Event> MultiBlockAnalysisFn::GetAllEvents() {

  vector<I3Event> allEvents;
  vector<I3Event> tempVect;

  for (int i=0; i<int(analysisFnVect_.size()); ++i) {
    analysisFnVect_[i]->PrepareAnalysis();
    tempVect.clear();
    tempVect = analysisFnVect_[i]->GetEventVector();

    for (int j=0; j<int(analysisFnVect_[i].size()); j++) {

      allEvents.push_back(tempVect[j]);
    }

  }

  return allEvents;

}

// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

void MultiBlockAnalysisFn::GetFlareGuess(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms) {


  // Note: already have Coord * srcCoord_ set;

  vector<I3Event> evs = GetAllEvents();

  vector<double> tVectortemp;
  vector<double> tVectorclose;

  double llhMax= -100.;
  double llhtemp, rms, avgtime;

  double sProb, eMaxRatio;
  const EnergyProb* eProb(NULL);
  I3Event event;

  for (int j=0; j<int(evs.size()); j++) {

    event = evs[j];

    sProb = event->ProbFrom(*srcCoord_);
    bProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event); //HMMM
    eProb = event->GetEnergyProbFn();
    eMaxRatio = eProb->GetEnergyMaxRatio(*event);


    if ( (eMaxRatio * sProb/bProb) > 10. ) {
        tVectorclose.push_back(event.GetMJD());
    }
  }

  sort(tVectorclose.begin(),tVectorclose.end());

  // I did some checks to see if you might find more flares if we tested a soft
  // spectrum for the seed as well, it didn't seem to be that useful.

//  const int ngamma=2;
//  double g[ngamma] = {2.0, 3.5};
  const int ngamma=1;
  double g[ngamma] = {2.0};

  int kmax = (int) close_;
  if (8 < kmax) { kmax = 8; }

  int n[] =  {1,2,3,4,5,10,15,20,25};
  int ns[] = {1,2,3,4,5, 5, 5, 5, 5};
  double p[4];

  int kmin=0;

  while ( tVectorclose.size()*1.0 < n[kmax] ) { kmax--; }

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

//      if (useE) {
        for (int h=0; h<ngamma; h++) {
          p[0] = ns[k]+1.;
          p[1] = g[h];
          p[2] = avgtime;
          p[3] = log10(rms);
          llhtemp = EvaluateLlh( p );
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_mean = avgtime;
            Guess_rms = rms;
            Guess_nsrc = ns[k]+1.;
          }
        }
        tVectortemp.clear();
//        } else { // I have no idea why you'd want to do the analysis without energy, but...
//         llhtemp = EvaluateLlh( ns[k]+1., 0., avgtime, log10(rms));
//         if (llhtemp > llhMax) {
//           llhMax = llhtemp;
//           Guess_mean = avgtime;
//           Guess_rms = rms;
//           Guess_nsrc = ns[k]+1.;
//         }
//         tVectortemp.clear();
//       }
    }
  }

  tVectorclose.clear();

  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, log10(Guess_rms));
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }
  } */

// Here are two quick functions for doing a lightcurve-based timedep
// analysis on multiple years/datasets. The initial conditions need to be
// calculated with all the data, the first one calculates a best possible
// lag btw photons and neutrinos and the second gets the max height of the
// lc and initial seed.


double MultiBlockAnalysisFn::SearchForLag(double laglimit) {

  double llhMax=-100.;
  double lagb=0., llhtemp;

  double p[] = {2, 2., 0., 0.};

  for (double d=-1.0*laglimit; d<laglimit; d=d+laglimit/10.){
    p[2] = d;
    llhtemp = EvaluateLlh( p );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      lagb = d;
    }
  }

  return lagb;
}

void MultiBlockAnalysisFn::SearchBlockSpace(string blocksFile, double lag, double & initV, double & maxT) {

  //string BlocksTimeFile_ = analysisFnVect_[0]->GetBlocksFile();

  //maxT = (GetHighestBlock(blocksFile.c_str()) + GetSecondHighestBlock(blocksFile_.c_str()))/2.;

  ifstream fin;
  fin.open(blocksFile.c_str());
  double blockbegin, blockdur, blocklev;
  double max=-1;
  double secondmax=-2;
  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;
    if (blocklev > secondmax) {
      if (blocklev > max) {
        secondmax = max;
        max = blocklev;
      } else {
        secondmax = blocklev;
      }
    }
  }
  fin.close();

  maxT = (max + secondmax) / 2.;

  double step = maxT/20.;
  double llhMax=-100.;
  double llhtemp;

  double p[] = {2, 2., lag, 0.};

  for (double d=0.;d<maxT;d+=step) {
    p[3] = d;
    llhtemp = EvaluateLlh( p );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      initV = d;
    }
  }

} //*/
