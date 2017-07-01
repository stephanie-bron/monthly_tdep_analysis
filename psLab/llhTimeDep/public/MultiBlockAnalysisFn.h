#ifndef LLH_MULTIBLOCKANALYSISFN_H_
#define LLH_MULTIBLOCKANALYSISFN_H_

//#include <exception>
#include "TFitterMinuit.h"

#include "rootExt/public/generalfunctions.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llhTimeDep/public/BlockLevel.h"
#include "llhTimeDep/public/NewLlhBlockTime.h"

class FluxBase;
class MultiBlockAnalysisFCN;


class MultiBlockAnalysisFn : public AnalysisFn {
 protected:
  TFitterMinuit* minuit_;
  MultiBlockAnalysisFCN* fcn_;

  vector<AnalysisFn*> analysisFnVect_;
  const ParTranslator* parTrans_;
  vector<MinuitParDef> parDefVect_;

  double logLambdaBest_;
  void StoreLogLambdaBest();

 public:

  int nPar;
  double laglimit;
  string blocksFile;
  int Ndof;

  MultiBlockAnalysisFn();
  virtual ~MultiBlockAnalysisFn() {
    if (minuit_) delete minuit_;
  }

  TFitterMinuit* Minuit() { return minuit_; }

  virtual void AddAnalysisFn(NewLlhBlockTime* llh) {

    for (int j=0;j <= (int) analysisFnVect_.size(); j++) {
        if (llh->ndof_!=Ndof && Ndof!=-1){
            cout << "ERROR: the newly added llh has a different number of ndof !!! stopping." << endl;
            exit(0);
        }
    }
    Ndof=llh->ndof_;
    cout << "for calculating p values will use Ndof= " << Ndof << endl;
    analysisFnVect_.push_back((AnalysisFn*) llh);
  }

  virtual void SetAnalysisSet(AnalysisSet*) {
    log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");
  }

  virtual void SetSearchCoord(const Coord& coord) {
    srcCoord_ = &coord;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->SetSearchCoord(coord);
    }
  }

  virtual void SetParDefs(vector<MinuitParDef>& parDefVect) {
    parDefVect_ = parDefVect;
  }

  virtual void SetParTranslator(const ParTranslator* pt) { parTrans_ = pt; }


  virtual void PrepareAnalysis() {
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->PrepareAnalysis();
    }
  }

  virtual void MaximizeLlh() {

    // Wipes out parameters (minuit crashes easily when changing existing pars)
    minuit_->Clear();
    minuit_->CreateMinimizer();  // default is kMigrad

    PrepareAnalysis();

    parDefVect_.clear();
    parDefVect_.push_back( MinuitParDef("nSrc",2.5,1., 0.,100.) );
    parDefVect_.push_back( MinuitParDef("gamma",2.5,0.1, 1., 4.) );

    double lagGuess = SearchForLag(laglimit);
    parDefVect_.push_back( MinuitParDef("lag", lagGuess, laglimit/10., -1.0*laglimit, laglimit) );

    double threshMax = ( GetHighestBlock(blocksFile) +
                         GetSecondHighestBlock(blocksFile) ) /2.;

    double guessThresh;
    SearchBlockSpace(blocksFile, lagGuess, guessThresh, threshMax);
    parDefVect_.push_back( MinuitParDef("thresh",guessThresh, threshMax/20., 0., threshMax) );

    for (int i=0; i<int(parDefVect_.size()); ++i) {

      const MinuitParDef& pd = parDefVect_[i];
      minuit_->SetParameter(i, pd.name.c_str(), pd.initValue, pd.initStepSize,
			    pd.lowLimit, pd.upLimit);
      //cout << pd.name.c_str() << " Init: " << pd.initValue << " max: " << pd.upLimit << endl;
    }
    minuit_->Minimize();
    StoreLogLambdaBest();
  }

  virtual double EvalFCN(const vector<double>& parVect) const;

  virtual double EvaluateLlh(double *parValueArray);

  virtual double GetPar(int i) const { return minuit_->GetParameter(i); }
  virtual const char* GetParName(int i) const { return minuit_->GetParName(i); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
  virtual double GetEstProb() const {
    double chiSq = 2. * Get_logLambdaBest();
    double p_temp, p;
    chisq_prob(chiSq, Ndof, &p_temp, &p);
    return p / 2.;  // one-sided chi-sq prob
  }

  //  vector<I3Event> GetAllEvents();
  //  void GetFlareGuess(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms);

  double SearchForLag(double laglimit);
  void SearchBlockSpace(string blocksFile, double lag, double & initV, double & maxT);

};



class MultiBlockAnalysisFCN : public ROOT::Minuit2::FCNBase {
 private:
  MultiBlockAnalysisFn* ptr;
 public:
  // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
  // This is what gets minimized: you have to define your likelihood
  virtual double operator() (const vector<double>& par) const;

  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<MultiBlockAnalysisFn*>(fn);
  }
};


#endif // LLH_MULTIBLOCKANALYSISFN_H_

