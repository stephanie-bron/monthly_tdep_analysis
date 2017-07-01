#include "llhTimeDep/public/MultiGaussAnalysisFn.h"
#include "rootExt/public/FunctionsRoot.h"

double MultiGaussAnalysisFCN::operator() (const vector<double>& parVect) const {
  double result = 0.;
  double rr;
  //cout << "parVect " << parVect[0] << " " << parVect[1] << " " << parVect[2] << " " << parVect[3] << endl;
  for (int i=0; i<int(ptr->analysisFnVect_.size()); ++i) {
    vector<double> individualParVect = ptr->parTrans_->Translate(i, parVect);
    rr=ptr->analysisFnVect_[i]->EvalFCN(individualParVect);
    result += rr;
    //Now we correct for the marginalization. 
    //Since it was applied N times on each individual llh
    //we have to remove it (N-1) times. 
    //Since result = -LogLikelihood it needs to be added (+). 
    //All marginalization weights are identical
  

    if (i>0 && ptr->analysisFnVect_[i]->JimsTerm_ && ptr->analysisFnVect_[i]->Get_margWeight() < 0.) { 
        result += ptr->analysisFnVect_[i]->Get_margWeight();
    }
  }
  return result;
}

MultiGaussAnalysisFn::MultiGaussAnalysisFn() {
    fcn_ = new MultiGaussAnalysisFCN();
    fcn_->Point(this);
    
    nPar = 4;
    histoForProb_ = false;
    seedWtMin=1000;
    
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

void MultiGaussAnalysisFn::StoreLogLambdaBest() {
    Double_t amin, edm, errdef;
    Int_t nvpar, nparx;
    minuit_->GetStats(amin, edm, errdef, nvpar, nparx);
    logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

    if ( analysisFnVect_[0]->monitorLevel_ > 1 ) { 
        printf("LogLambdaBest=%12.6lg\n",logLambdaBest_ );
    } 

    // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
    // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
    // negative. Since this will cause probability calculation to choke, we 
    // fix it here
  
    if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}

vector<I3Event> MultiGaussAnalysisFn::GetAllEvents() {
    vector<I3Event> allEvents;
    vector<I3Event> tempVect;

    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
        NewLlhGausTime* i3an = dynamic_cast<NewLlhGausTime*>(analysisFnVect_[i]);
    
        if (i3an->Get_nEvents()==0) { i3an->PrepareAnalysis(); }
        tempVect.clear();
        tempVect = i3an->GetEventVector();
    
        for (int j=0; j<int(tempVect.size()); j++) {
        allEvents.push_back(tempVect[j]);
        }
    }
    return allEvents;
}


void MultiGaussAnalysisFn::MaximizeLlh() {
    // Wipes out parameters (minuit crashes easily when changing existing pars)
    minuit_->Clear();
    minuit_->CreateMinimizer();  // default is kMigrad

    PrepareAnalysis();
    parDefVect_.clear();

    //find the bounds on spectral index  
    gammaMin_ = analysisFnVect_[0]->GetGammaMin();
    gammaMax_ = analysisFnVect_[0]->GetGammaMax();
    
    gammaMin_ = 1.;
    gammaMax_ = 4.;
    
    GetFlareGuessGauss(nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_, sigmamin_);
    if (gammaGuess_<1.) gammaGuess_=1.1;
    if (gammaGuess_>4.) gammaGuess_=3.9;
    if (sigmamin_<1E-8) sigmamin_=1E-8;
    double nsrcMax=nEventsTot_;
    
    for (unsigned int i=0; i<int(analysisFnVect_.size()); ++i) {
        double maxw=-1;
        for (int k=0;k<(parTrans_->srcWeightVect_[i]).size();k++) {
            if (maxw<parTrans_->srcWeightVect_[i][k]) maxw=parTrans_->srcWeightVect_[i][k];
        }
        //cout << "analysisFnVect_[i]->nEventsTot_ " << analysisFnVect_[i]->nEventsTot_ << " " << maxw << " " << analysisFnVect_[i]->nEventsTot_/maxw << " "<< nsrcMax << endl;
        nsrcMax= min(analysisFnVect_[i]->nEventsTot_/maxw,nsrcMax);
    }    
    
    
    parDefVect_.push_back( MinuitParDef("nSrc",nsrcGuess_,0.1, 0.,0.5*nsrcMax) );
    parDefVect_.push_back( MinuitParDef("gamma",gammaGuess_,0.1, gammaMin_, gammaMax_) );   
    parDefVect_.push_back( MinuitParDef("mean", meanGuess_, sigmaGuess_, tmin_, tmax_) );

    // The width is fit in log10 space, following the untriggered flare search.
    // It probably doesn't make too much of a difference here, since it only searches over a few
    // orders of magnitude in width.
    
    parDefVect_.push_back( MinuitParDef("sigma",log10(sigmaGuess_), 1., log10(sigmamin_), log10((tmax_-tmin_)/2.))  );
    //cout << "par set ss " << sigmaGuess_ << " " << sigmamin_ << " " << tmax_ << " " << tmin_ << endl;
    //parDefVect_.push_back( MinuitParDef("nSrc", 6, 0.1, 0, 169301));
    //parDefVect_.push_back( MinuitParDef("gamma", 2.4, 0.1, 1, 4));
    //parDefVect_.push_back( MinuitParDef("mean", 56608.7, 111.784, 56062.4, 57160));
    //parDefVect_.push_back( MinuitParDef("sigma", 2.04838, 1, -3.63553, 2.72327));

        
    for (int i=0; i<int(parDefVect_.size()); ++i) {
        const MinuitParDef& pd = parDefVect_[i];
        minuit_->SetParameter(i, pd.name.c_str(), pd.initValue, pd.initStepSize,pd.lowLimit, pd.upLimit);
        //cout << "par set " << i<< " " << pd.name.c_str()<< " " << pd.initValue<< " " << pd.initStepSize<< " " <<pd.lowLimit<< " " << pd.upLimit << endl;
    }
    minuit_->Minimize();
    StoreLogLambdaBest();
    
    nSrcBest_  = minuit_->GetParameter(0);
    gammaBest_ = minuit_->GetParameter(1);
    meanBest_  = minuit_->GetParameter(2); 
    sigmaBest_ = pow( 10., minuit_->GetParameter(3) );
}

void MultiGaussAnalysisFn::GetFlareGuessGauss(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms, double & sigmamin_){

    vector<double> tVectortemp;
    vector<double> tVectorclose;
    vector<I3Event> evs;
    vector<double> LCBkgProb; //first it will be set to 1. if no LCBkgProb to be used and else -1 for normal LCBkgProb, -2 for folded
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
        vector<I3Event> eVect=analysisFnVect_[i]->GetEventVector();
        evs.insert(evs.end(), eVect.begin(), eVect.end());
        for (int ii=0;ii<eVect.size();ii++){
            if (analysisFnVect_[i]->useLCBkgProb_) {
                if ( analysisFnVect_[i]->UseFolded_ ) LCBkgProb.push_back(analysisFnVect_[i]->lcBkgProb_->BackgroundLCProb_folded(eVect[ii]));
                else LCBkgProb.push_back(analysisFnVect_[i]->lcBkgProb_->BackgroundLCProb(eVect[ii]));
            }                                                  
            else LCBkgProb.push_back(1.);
        }
    }
    
    vector<double> ttt;

    for (unsigned int k=0;k<evs.size();k++) ttt.push_back( evs[k].GetMJD());
    sort(ttt.begin(),ttt.end());
    double timed;
    sigmamin_=1000.;
    for (unsigned int i=1; i<ttt.size(); i++) {//This assumes events are time-ordered.
        timed = ttt[i] - ttt[i-1];
        if (analysisFnVect_[0]->monitorLevel_ > 2)  cout << ttt[i] << " "; 
        if (fabs(timed) < sigmamin_) sigmamin_ = fabs(timed);
    }
    if (analysisFnVect_[0]->monitorLevel_ > 2) cout << endl;
    sigmamin_ /= sqrt(2.);
    if (analysisFnVect_[0]->monitorLevel_ > 2) { cout << " sigmamin= " << sigmamin_ << endl; }
    ttt.clear();
    
    
    double llhtemp, rms, avgtime;
    double sProb,bProb,eMaxRatio;
    double llhMax= -100.;

    for (int j=0; j<int(evs.size()); j++) { 
        
        sProb = evs[j].ProbFrom(*srcCoord_);
        bProb = evs[j].GetBkgSpaceProbFn()->GetBkgProbDensity(evs[j]);
        eMaxRatio = evs[j].GetEnergyProbFn()->GetEnergyMaxRatio(evs[j]);
        bProb = bProb*LCBkgProb[j];
        if ( (eMaxRatio * sProb/bProb) > seedWtMin ) tVectorclose.push_back(evs[j].GetMJD());           
    }
    if (tVectorclose.size() <10) {
        tVectorclose.clear();
        for (int j=0; j<int(evs.size()); j++) tVectorclose.push_back(evs[j].GetMJD());
    }

    sort(tVectorclose.begin(),tVectorclose.end());

    if (analysisFnVect_[0]->monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". Sigmamin is " << sigmamin_ << endl; }
  
    const int ngamma=1;
    double g[ngamma] = {2.0};

    int kmax = 10; // MUST BE <=11

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
            for (unsigned int ii=0;ii<tVectortemp.size();ii++) sum += tVectortemp[ii];
            avgtime = sum/tVectortemp.size();

            for (unsigned int ii=0;ii<tVectortemp.size();ii++) sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
            //calculate the rms to use as the sigma
            rms = sqrt( sumsq/tVectortemp.size() );

            //cout << " GUESS k " << k << " i " << i << " avgtime " << avgtime << " rms " << rms << endl;

            llhtemp = 0.;

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
        }
    }

    double sllhmax = -100;
    double sllhtemp;
    for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
        sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, Guess_rms );
        if (sllhtemp > sllhmax ) {
            Guess_gamma = d;
            sllhmax = sllhtemp;
        }
    }
    //cout << " GUESS " << Guess_nsrc<< " " << Guess_gamma<< " " << Guess_mean<< " " << Guess_rms<< " " << sigmamin_ << endl;
}
    
double MultiGaussAnalysisFn::EvaluateLlh(double *parValueArray) {
    vector<double> parVect;
    for (int i=0; i<nPar; ++i) {
        parVect.push_back(parValueArray[i]);
    }
    double minusLlh = EvalFCN(parVect);
    return -minusLlh;   // that is, max llh = - (minimizer result)
}

double MultiGaussAnalysisFn::EvalFCN(const vector<double>& parVect) const {

  return (*fcn_)(parVect);
}

double MultiGaussAnalysisFn::GetProbFromHisto(double teststat) const { 
    int bin = pvalHisto_->FindBin(teststat);
    double ptemp = pvalHisto_->GetBinContent(bin);
    return ptemp;
}

void MultiGaussAnalysisFn::SetNullTestStat(TH1D * inputhisto) {
  
  // This reads in a specific TH1D as the null test statistic
  // to use for p-values instead of using a chisquare distribution.
  // It also used to fit an exponential to the upper tail (default set to top 0.01%).
  
  // It takes the raw test statistic distribution (pdf) as input
  // and then makes the cdf to work with later.
  
    histoForProb_ = true;
  
    pvalHisto_ = new TH1D();
    nullTestStat_ = new TH1D();
  
    char sc[] = "Scale";
    pvalHisto_ = DescendingCumulate(inputhisto,sc);
    nullTestStat_ = (TH1*)inputhisto->Clone("hnew");  
}