
#include "llh/public/MultiAnalysisFn.h"
#include "iostream"



double MultiAnalysisFCN::operator() (const double *par) const {
  double result = 0.;
  

  for (int i = 0; i < int(ptr->analysisFnVect_.size()); ++i) {

    for (int j = 0; j <  ptr->nPar_; j++)
      {
	if (par[j] != par[j]){
	  log_warn("MultiAnalysisFn: Global parameter is NaN. Returning 1e50 in the minimizer");
	  
	  return 1e50;
	}
      }
    

    //Create a vector<double> it is needed for the Translator
    vector<double> parVect;
    for (int j = 0; j <  ptr->nPar_; j++)
      {
	parVect.push_back(par[j]);
      }
    
    vector<double> individualParVect = ptr->parTrans_->Translate(i, parVect);
    
    //The individual llh still use EvalFCN(vector<double> par) instead of EvalFCN(double *par)
    result += ptr->analysisFnVect_[i]->EvalFCN(individualParVect);
    
  }
  
  return result;
}



MultiAnalysisFn::MultiAnalysisFn() {

  useSimplex_ = false;
  
  fcn_ = new MultiAnalysisFCN();
  fcn_->Point(this);
  
  //Default is migrad
  migrad_  = ROOT::Math::Factory::CreateMinimizer("Minuit2", "MIGRAD");
  simplex_  = ROOT::Math::Factory::CreateMinimizer("Minuit2", "SIMPLEX");
  
}


double MultiAnalysisFn::EvalFCN(const double *parVect) const {
  return (*fcn_)(parVect);
}



double MultiAnalysisFn::EvaluateLlh(double *parValueArray) {

  double minusLlh = EvalFCN(parValueArray);

  return -minusLlh;   // that is, max llh = - (minimizer result)
  
}


void MultiAnalysisFn::StoreLogLambdaBest() 
{

  if(migrad_->Status() == 1 && useSimplex_)
    logLambdaBest_ = -simplex_->MinValue();
  else
    logLambdaBest_ = -migrad_->MinValue();
  
  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}

