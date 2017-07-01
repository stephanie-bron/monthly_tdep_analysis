#ifndef LLH_MULTIGAUSSANALYSISFN_H_
#define LLH_MULTIGAUSSANALYSISFN_H_

#include "llh/public/MinuitAnalysisFn.h"
#include "llhTimeDep/public/NewLlhGausTime.h"
#include "llh/public/LlhFunctionsBase.h"
#include "rootExt/public/generalfunctions.h"

//#include "llh/public/I3Event.h"

#include "TFitterMinuit.h"

class FluxBase;
class MultiGaussAnalysisFCN;

class MultiGaussAnalysisFn : public AnalysisFn {
    protected:
        TFitterMinuit* minuit_;
        MultiGaussAnalysisFCN* fcn_; 

        vector<NewLlhGausTime*> analysisFnVect_; // this needs to be specifically
                                              // for periodic analyses, not an <I3Analysis>
        const NewLlhGausTime_ParTranslator* parTrans_;
        vector<MinuitParDef> parDefVect_;
        double tmin_;
        double tmax_;

        double nEventsTot_;

        double nsrcGuess_;
        double gammaGuess_;
        double meanGuess_;
        double sigmaGuess_;
        double sigmamin_;

        double gammaMin_;
        double gammaMax_;
                
        double logLambdaBest_;
        void StoreLogLambdaBest();

        double nSrcBest_;
        double gammaBest_;
        double meanBest_;
        double sigmaBest_;

        bool histoForProb_;
        TH1 * pvalHisto_;   // distribution to estimate p-values
        TH1 * nullTestStat_; // a few parameters for using a custom teststat

    public:
        double seedWtMin;
        int nPar;

        MultiGaussAnalysisFn();
        virtual ~MultiGaussAnalysisFn() { 
            if (minuit_) delete minuit_;
        }

        TFitterMinuit* Minuit() { return minuit_; }
        
        virtual void AddAnalysisFn(AnalysisFn* llh) {
            NewLlhGausTime* llh1 = dynamic_cast<NewLlhGausTime*>(llh);
            analysisFnVect_.push_back(llh1);
        }
        
        virtual void SetAnalysisSet(AnalysisSet*) {log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");}

        virtual void SetSearchCoord(const Coord& coord) {
            srcCoord_ = &coord;
            for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->SetSearchCoord(coord);
            }
        }
        
        virtual void GetFlareGuessGauss(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms, double & sigmamin_);
        
        virtual void SetParDefs(vector<MinuitParDef>& parDefVect) { parDefVect_ = parDefVect;}
        virtual void AddParDef(MinuitParDef parDef){
            parDefVect_.push_back(parDef);
            }
        virtual void SetParTranslator(const NewLlhGausTime_ParTranslator* pt) { parTrans_ = pt; }

        virtual void PrepareAnalysis() {
            nEventsTot_ = 0.;
            for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->PrepareAnalysis();
                nEventsTot_ += analysisFnVect_[i]->Get_nEvents();
            }
        }
                
        virtual void MaximizeLlh();
        
        vector<I3Event> GetAllEvents();
        virtual double EvaluateLlh(double *parValueArray);
        double EvaluateLlh( double a, double b, double c, double d) {    
            double dd[] = {a, b, c, log10(d) }; //use log10(sigma) in minimizer
            double Llh = EvaluateLlh(dd);
            return Llh;
        }
        double Get_nSrcBest()  { return nSrcBest_;}
        double Get_gammaBest() { return gammaBest_;}
        double Get_meanBest()  { return meanBest_;}
        double Get_sigmaBest() { return sigmaBest_;}

        double GetNsrcGuess() { return nsrcGuess_; }
        double GetMeanGuess() { return meanGuess_; }
        double GetSigmaGuess(){ return sigmaGuess_;}

        virtual double EvalFCN(const vector<double>& parVect) const;
        
        double GetProbFromHisto(double teststat) const;
        void SetNullTestStat(TH1D * inputhisto);

        virtual double GetPar(int i) const { return minuit_->GetParameter(i); }
        virtual double Get_logLambdaBest() const { return logLambdaBest_; }
        virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
        virtual double GetSigmaMin() { return sigmamin_;}
        virtual double GetEstProb() const {
            if (histoForProb_)  return GetProbFromHisto( Get_logLambdaBest() ); 
            double chiSq = 2. * Get_logLambdaBest();
            if(chiSq<0)  chiSq=0.; 
            double p_temp, p;
            int nDoF = analysisFnVect_[0]->ndof_;
            chisq_prob(chiSq, nDoF, &p_temp, &p);
            return p / 2.;  // one-sided chi-sq prob
        }
  
        void SetTimeBounds(double tmin, double tmax) {
            tmin_ = tmin;
            tmax_ = tmax;
        }
};

class MultiGaussAnalysisFCN : public ROOT::Minuit2::FCNBase {
    private:
        MultiGaussAnalysisFn* ptr;
    public:
        // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
        // This is related to how "errors" are defined, see Minuit2 for documentation
        virtual double Up() const { return 0.5; }

        // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
        // This is what gets minimized: you have to define your likelihood
        virtual double operator() (const vector<double>& par) const;

        // Here's where the connection is made to FCN can access the data
        virtual void Point(AnalysisFn* fn) {
            ptr = dynamic_cast<MultiGaussAnalysisFn*>(fn); 
        }
};


#endif // LLH_MultiGAUSSANALYSISFN_H_
