#ifndef LLH_I3SIGNALGENERATOR_H_
#define LLH_I3SIGNALGENERATOR_H_

#include <vector>

#include "rootExt/public/log_report.h"

//#include "fluxus/public/FluxFunction.h"
// Forward Declarations (when feasible, more efficient than including headers)
class FluxBase;

#include "llh/public/I3Event.h"
#include "llh/public/TimePdf.h"


//
// ABSTRACT BASE CLASS
//



class I3SignalGenerator : public SourceModule {
 public:
  virtual ~I3SignalGenerator() { }

  virtual I3SignalGenerator* Clone() const = 0;
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  // "SET" FUNCTIONS

  virtual void SetLivetime(double livetime) = 0;

  // "GET" FUNCTIONS

  virtual double GetLivetime() const = 0;

  virtual double GetMeanSrcNev() const = 0;
  // mean number of signal events generated for flux model "installed"

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const = 0;
  // mean number of signal events... for any other model specified

  // how much to scale model flux up or down to match Nev
  virtual double GetFluxScaleForNev(double Nev) const {
    double meanNev = GetMeanSrcNev();
    if (meanNev>0) { return Nev / meanNev; }
    return 0.;
  }

  virtual void SetTimePdf(TimePdf * timePdf) {
    if ( timePdf->GetLivetime() ) { }
    cout << "setTimePdf not implemented!\n";
  }
  virtual TimePdf * GetTimePdf() const = 0;

  // GENERATOR

  virtual int GenerateEventEntry() = 0;
  // Just get *which* event would be generated, not the event itself
  // (useful for diagnosing the statistics of the signal simulation files)

  virtual I3Event GenerateEvent() = 0;
};



//
// MULTI-SIGNAL GENERATOR: 
//
// Combine different I3SignalGenerator objects into one signal.
//
//
// You can add as many I3SignalGenerator's of any type as you like.
// When you GenerateEvent() with Multi, it will randomly choose between 
// the different sources, based on their relative weighting.

// To begin with, their weighting is based on their event normalizations: 
// GetMeanSrcNev().  This depends only on how the source itself was 
// configured.  There are two ways this can be altered in Multi:

// 1) change the enhanceFactor for the source in Multi, which will 
// multiply the GetMeanSrcNev normalization.

// 2) the more obscure way is to change the livetime of one of the sources.  
// Each source is added with the livetime it started with, but you can 
// change it later with AccessSignalPtr(i)->SetLivetime(livetime).  (You
// can in fact get complete access to the source with this method.)

// Note that the function I3MultiSignalGenerator::SetLivetime(livetime) will 
// reset all signals with the same global livetime.


class I3MultiSignalGenerator : public I3SignalGenerator {
 private:
  double livetime_;
  vector<I3SignalGenerator*> signalPtrVect_;
  vector<double> enhanceFactorVect_;
  TimePdf * timePdf_; // with blocks and such needing to be smart about 
                      // normalizations, it sounds like a good idea to keep
                      // an extra copy here.

 public:
  I3MultiSignalGenerator() : livetime_(0.), timePdf_(NULL) { }

  virtual ~I3MultiSignalGenerator() {
    for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      delete signalPtrVect_[i];
    }
  }

  // Need a deep copy, since new copies must be made for signalPtr's
  virtual I3SignalGenerator* Clone() const;

  // "SET" FUNCTIONS

  // This sets all livetimes simultaneously;
  // otherwise they can be independent of each other (e.g. a burst?)
  virtual void SetLivetime(double livetime);
  
  void SetTimePdf(TimePdf * tPdf);

  void AddSignal(const I3SignalGenerator& signal, 
		    double enhanceFactor = 1.);

  void SetEnhanceFactor(unsigned int i, double factor) {
    if (i<enhanceFactorVect_.size()) { enhanceFactorVect_[i] = factor; }
    else { log_error("%d enhanceFactor not assigned.\n",i); }
  }

  // "GET" FUNCTIONS

  virtual double GetLivetime() const {return livetime_;}
  
  TimePdf * GetTimePdf() const {return timePdf_;}

  virtual double GetMeanSrcNev() const;
  // adds up enhanceFactorVect_[i] * signalVect_[i].GetMeanSrcNev()

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;
    // same as above, for any model
  
  virtual double GetMeanSrcNevBase() const;
  // adds up signalVect_[i].GetMeanSrcNev()

  virtual double GetMeanSrcNevForFluxModelBase(const FluxBase& fluxModel) const;
    // same as above, for any model


  I3SignalGenerator* AccessSignalPtr(unsigned int i) {
    if (i<signalPtrVect_.size()) { return signalPtrVect_[i]; }
    else { log_warn("%d signalPtr not assigned.\n",i); return NULL; }
  }

  double GetEnhanceFactor(unsigned int i) {
    if (i<enhanceFactorVect_.size()) { return enhanceFactorVect_[i]; }
    else { log_warn("%d enhanceFactor not assigned.\n",i); return 0.; }
  }

  // GENERATOR

  int GenerateEventEntry() {
    log_fatal("Entry numbers are not well defined and currently not"
	      "implemented in I3MultiSignalGenerator.\n");
    return -1;  // what else can you do?
  }

  I3Event GenerateEvent();
};
  


// 
// SIMPLE GENERATOR FOR A POINT SOURCE
//



class I3PointGenerator : public I3SignalGenerator {
 private:
  double livetime_;
  double weightSum_;
  double weightMax_;
  EquatorialDeg sourceCoord_;
  vector<I3Event> candidateEvents_;

  int timeAzBins_; // to make sure we're injecting signal
                   // with the right local coordinate structure
    
  TimePdf * sourceTimePdf_;
  
  double fluenceSum_;

 public:

  I3PointGenerator();
  I3PointGenerator(const vector<I3Event>& inputEvents, 
		   const FluxBase& fluxModel,
		   const EquatorialDeg& sourceCoord, 
		   double livetime);
  I3PointGenerator(const vector<I3Event>& inputEvents, 
		   const FluxBase& fluxModel,
		   const EquatorialDeg& sourceCoord, 
		   double livetime,
		   TimePdf * timepdf);
  ~I3PointGenerator() { }

  I3SignalGenerator* Clone() const { return new I3PointGenerator(*this); }
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class, thusly:
  // { return new DerivedClass(*this); }

  // "SET" FUNCTIONS

  // If you set events, you have to set weights too...
  void SetLivetime(double livetime) { livetime_ = livetime; }

  void SetCandidateEvents(const vector<I3Event>& inputEvents, 
			  const FluxBase& fluxModel,
			  const EquatorialDeg& sourceCoord);

  // ... but you can change flux or livetime without changing anything else
  void SetCandidateWeights(const FluxBase& fluxModel);

  void SetTimePdf(TimePdf * tPdf) {
    sourceTimePdf_ = tPdf->Clone();
    sourceTimePdf_->fillHisto(1e3);
  }
  void SetTimeAzBins(int nbins) { timeAzBins_ = nbins; }
  
  // "GET" FUNCTIONS

  
  
  virtual double GetLivetime() const {return livetime_;}
  virtual double GetMeanSrcNev() const {return weightSum_ * livetime_;}
  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;
  TimePdf * GetTimePdf() const {return sourceTimePdf_;}
  double GetFluenceSum() {return fluenceSum_;}
  
  const vector<I3Event>* GetCandidateEvents() const {return &candidateEvents_;}
  EquatorialDeg GetSourceCoord() const {return sourceCoord_;}
  double GetWeightSum() const {return weightSum_;}
  double GetWeightMax() const {return weightMax_;}

  // GENERATOR

  int GenerateEventEntry();
  I3Event GenerateEvent();
};


#endif // LLH_I3SIGNALGENERATOR_H_
