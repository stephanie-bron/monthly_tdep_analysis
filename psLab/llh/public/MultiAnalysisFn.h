#ifndef LLH_MULTIANALYSISFN_H_
#define LLH_MULTIANALYSISFN_H_


#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"



#include "rootExt/public/generalfunctions.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

class FluxBase;
class MultiAnalysisFCN;


class MultiAnalysisFn : public AnalysisFn {
 protected:


  ROOT::Math::Minimizer* migrad_;
  ROOT::Math::Minimizer* simplex_;

  MultiAnalysisFCN* fcn_; 

  vector<AnalysisFn*> analysisFnVect_;

  const ParTranslator* parTrans_;
  vector<MinuitParDef> parDefVect_;
  
  int status_;
  
  double logLambdaBest_;
  void StoreLogLambdaBest();

 public:
  
  bool useSimplex_;
  
  int nPar_;

  void SetUseSimplex(bool useSimplex)
  {
    useSimplex_ = useSimplex;
  };

  MultiAnalysisFn();
  virtual ~MultiAnalysisFn() { 
    if (migrad_) delete migrad_;
    if (simplex_) delete simplex_;
  }
    
  virtual void AddAnalysisFn(AnalysisFn* llh) {
    analysisFnVect_.push_back(llh);
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
    nPar_ = parDefVect_.size();
  }

  virtual void SetParTranslator(const ParTranslator* pt) { parTrans_ = pt; }


  virtual void PrepareAnalysis() {
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->PrepareAnalysis();
    }
  }

  virtual void MaximizeLlh() {
    
    PrepareAnalysis();
    
    
    if(nPar_ == 0)
      log_error("No parameters defined!");
    

    ROOT::Math::Functor f(*fcn_, nPar_);

    migrad_->SetFunction(f);
    simplex_->SetFunction(f);

    migrad_->Clear();  
    simplex_->Clear();
    
    migrad_->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
    simplex_->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
    
    // Call Migrad with 500 iterations maximum 
    migrad_->SetMaxIterations(500);  
    simplex_->SetMaxIterations(500);  

    migrad_->SetTolerance(0.001);
    simplex_->SetTolerance(0.001);

    
    for (int i=0; i<int(parDefVect_.size()); ++i) {
      
      const MinuitParDef& pd = parDefVect_[i];
      migrad_->SetLimitedVariable(i, pd.name.c_str(), pd.initValue, pd.initStepSize,
				  pd.lowLimit, pd.upLimit);
      simplex_->SetLimitedVariable(i, pd.name.c_str(), pd.initValue, pd.initStepSize,
			       pd.lowLimit, pd.upLimit);
    }
    
    // Set Print Level 
    // NOTE: this is overwritten by global variable gErrorIgnoreLevel
    // Please put gErrorIgnoreLevel = kWarning in your scripts otherwise you will get tons of INFO messages 
    migrad_->SetPrintLevel(-1);

    //Let's minimuze first only migrad
    migrad_->Minimize();

    
    //Change to SIMPLEX if migrad failed!
    status_ = migrad_->Status();
    
    if(status_ == 1 && useSimplex_) //Covar was made positive
      {
	simplex_->SetPrintLevel(0);
	simplex_->Minimize();
	status_ = simplex_->Status();
      }
    
    
    StoreLogLambdaBest();    
  }

  
  virtual double EvalFCN(const double *parVect) const;

  virtual double EvaluateLlh(double *parValueArray);

  virtual double GetPar(int i) const { 
    
    const double *xs;
    if(migrad_->Status() == 1 && useSimplex_)
      xs = simplex_->X();
    else
      xs = migrad_->X();

    return xs[i]; 
  }

  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
  virtual double GetEstProb() const {
    double chiSq = 2. * Get_logLambdaBest();
    double p_temp, p;
    int nDoF = migrad_->NFree();
    chisq_prob(chiSq, nDoF, &p_temp, &p);
    return p / 2.;  // one-sided chi-sq prob
  }

  int GetMinimizerStatus() { return status_; }
};



//A wrapper to the operator()
class MultiAnalysisFCN {

 private:
  MultiAnalysisFn* ptr;
 public:
  virtual double operator() (const double* par) const;
  
  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<MultiAnalysisFn*>(fn); 
  }
};



#endif // LLH_MULTIANALYSISFN_H_

