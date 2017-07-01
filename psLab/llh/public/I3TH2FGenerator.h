#ifndef I3TH2FGENERATOR_H_
#define I3TH2FGENERATOR_H_

#include <vector>
#include "rootExt/public/log_report.h"
#include "llh/public/I3SignalGenerator.h"
#include "llh/public/I3Event.h"
#include "llh/public/EventLoader.h"
#include "fluxus/public/FluxFunction.h"
#include "TH2F.h"

// Forward Declarations (when feasible, more efficient than including headers)

// Generator for a signal distributed according to a 2D map, given as a TH2F.

class I3TH2FGenerator : public I3SignalGenerator {
 private:
  double livetime_;
  double weightSum_;
  double weightMax_;
  vector<I3Event> candidateEvents_; 

  // NEW MEMBERS FOR G.P.
  EquatorialDeg testCoord_;

  // Keep pointers for on-the-fly candidate event selection
  EventLoader *evLoaderPtr_; // Not const since LoadSourceEvents is now not const
  const FluxBase *fluxModelPtr_;

  // Product of FluxModel and EmissionModel (yields signal event distribution)
  TH2F *hSrcWeightsTable_;

  // Determine range and binning of srcWeightsTable_
  double xMin_;
  double xMax_;
  int nXbins_;
  double yMin_;
  double yMax_;
  int nYbins_;

 public:
  I3TH2FGenerator();
  I3TH2FGenerator(EventLoader& evLoader,
		   const FluxBase& fluxModel,
		   double livetime);

  ~I3TH2FGenerator() { }

  I3SignalGenerator* Clone() const { return new I3TH2FGenerator(*this); }
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class, thusly:
  // { return new DerivedClass(*this); }

  // "SET" FUNCTIONS

  // If you set events, you have to set weights too...
  void SetLivetime(double livetime) { livetime_ = livetime; }

  // This now needs to be done for each event generated!!
  void SetCandidateEvents();

  void SetEventLoader(EventLoader& el) {
    evLoaderPtr_ = &el;
  }
			  
  void SetFluxModelPtr( const FluxBase& fluxModel ) {
    fluxModelPtr_ = &fluxModel;
  }
			  
  void SetCandidateWeights();

  // "GET" FUNCTIONS

  virtual double GetLivetime() const {return livetime_;}
  // JD: livetime already accounted for in construction...
  virtual double GetMeanSrcNev() const {return weightSum_;}
  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const{return weightSum_;};

  TimePdf * GetTimePdf() const {
    //cout << "Time Pdf not implemented!\n"; 
    return NULL;
  }

  const vector<I3Event>* GetCandidateEvents() const {return &candidateEvents_;}
  double GetWeightSum() {
	CalculateWeightSum();
	return weightSum_;
  }
  void SetWeightSum(double weightSum) {
	weightSum_ = weightSum;
  }
  double GetWeightMax() const {return weightMax_;}

  void CalculateWeightSum();

  EquatorialDeg GetTestCoord() const {return testCoord_;}

  void SetSrcWeightsTable(TH2F* hSrcWeights); 

  // Can set the binning by hand
  void SetBinning(double xMin, double xMax, int nXbins, double yMin, double yMax, int nYbins) {
    xMin_ = xMin;
    xMax_ = xMax;
    nXbins_ = nXbins;
    yMin_ = yMin;
    yMax_ = yMax;
    nYbins_ = nYbins;
  }
  // Can set the binning to match an existing histogram
  void SetBinning(TH2F* hBinning) { 
    printf("Using an existing histo to set binning\n");
    xMin_ = hBinning->GetXaxis()->GetXmin();
    xMax_ = hBinning->GetXaxis()->GetXmax();
    nXbins_ = hBinning->GetXaxis()->GetNbins();
    yMin_ = hBinning->GetYaxis()->GetXmin();
    yMax_ = hBinning->GetYaxis()->GetXmax();
    nYbins_ = hBinning->GetYaxis()->GetNbins();
  }

  // Pick a spot according to TH2F, stores as testCoord_
  void GenerateCoord();

  TH2F* GetSrcWeightsTable(); 

  // GENERATOR
  int GenerateEventEntry();
  I3Event GenerateEvent();
};


#endif // I3TH2FGENERATOR_H_
