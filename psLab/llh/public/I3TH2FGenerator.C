#include "rootExt/public/randomfunctions.h"
#include "llh/public/I3TH2FGenerator.h"
#include "fluxus/public/FluxFunction.h"
#include "llh/public/CoordinateTransform.h"
#include "TFile.h"

// I3TH2FGenerator

I3TH2FGenerator::I3TH2FGenerator() : 
  livetime_(0),
  weightSum_(0),
  weightMax_(0),
  xMin_(0),
  xMax_(360),
  nXbins_(720),
  yMin_(-10),
  yMax_(10),
  nYbins_(40)
{ }


I3TH2FGenerator::I3TH2FGenerator(
				   EventLoader& evLoader,
				   const FluxBase& fluxModel,
				   double livetime) :
  livetime_(livetime),
  weightSum_(0),
  weightMax_(0),
  xMin_(0),
  xMax_(360),
  nXbins_(720),
  yMin_(-10),
  yMax_(10),
  nYbins_(40)
{
  fluxModelPtr_ = &fluxModel;
  evLoaderPtr_ = &evLoader;
}

void I3TH2FGenerator::CalculateWeightSum() {
  double sum = 0.0;
  for (int xIndex=1; xIndex <= nXbins_; xIndex++) {
    for (int yIndex=1; yIndex <= nYbins_; yIndex++) {

      double RaDeg  = hSrcWeightsTable_->GetXaxis()->GetBinCenter(xIndex);
      double DecDeg = hSrcWeightsTable_->GetYaxis()->GetBinCenter(yIndex);
 
      double  enhancementFactor = hSrcWeightsTable_->GetBinContent(xIndex, yIndex);
      //VERY IMPORTANT!!!! these histograms are not normalized to sum to 1 but to per str per bin!!!
      enhancementFactor = enhancementFactor/hSrcWeightsTable_->Integral();
	
cout<<enhancementFactor<<endl; 
      if(enhancementFactor>0.0){ 
        candidateEvents_.clear();
       
        testCoord_ = EquatorialDeg(RaDeg, DecDeg);
        evLoaderPtr_->LoadSourceEvents(candidateEvents_, testCoord_);
       
        double binWeightSum = 0.0;
       
        for (vector<I3Event>::iterator e = candidateEvents_.begin(); e != candidateEvents_.end(); e++) {
          I3MCParameters mc = e->GetMCParams();
       
          double weight = fluxModelPtr_->GetFlux(mc.mcEnergy) * mc.PS_FlatSpectrumRate;
          // this gets multiplied by livetime for number of events, so
          // effectively it could be considered "number of events for 1 sec" here.
       
          mc.srcWeight = weight;
          e->SetMCParams(mc);
       
          binWeightSum += weight;
        } 
        
        sum += binWeightSum*enhancementFactor;
      }
    }
  }
  weightSum_ =  sum * livetime_;
}

// This now needs to be done for each event generated!!
void I3TH2FGenerator::SetCandidateEvents() {
  candidateEvents_.clear();
  GenerateCoord(); // overwrites testCoord_ with new random location
  evLoaderPtr_->LoadSourceEvents(candidateEvents_, testCoord_);
  SetCandidateWeights();
}


// Weights do not include the livetime (therefore, livetime can be changed 
// independently at any time).  Effectively, this means each weight stands
// for "the mean number of such events in 1 sec"
// (since sec is unit of livetime).

void I3TH2FGenerator::SetCandidateWeights() {
  weightMax_ = 0.;

  for (vector<I3Event>::iterator e = candidateEvents_.begin();
       e != candidateEvents_.end(); e++) {
    I3MCParameters mc = e->GetMCParams();

    double weight = fluxModelPtr_->GetFlux(mc.mcEnergy) * mc.PS_FlatSpectrumRate;
    // this gets multiplied by livetime for number of events, so
    // effectively it could be considered "number of events for 1 sec" here.

    if (weight<0) { 
      log_error("ERROR: Tried to assign weight<0\n"); 
      weight = 0.;
    }
      
    mc.srcWeight = weight;
    e->SetMCParams(mc);

    if (weight > weightMax_) weightMax_ = weight;
  }
}


void I3TH2FGenerator::SetSrcWeightsTable(TH2F *hSrcWeights) {

  assert(hSrcWeights->GetNbinsX() == nXbins_);
  assert(hSrcWeights->GetXaxis()->GetBinLowEdge(1) == xMin_);
  assert(hSrcWeights->GetXaxis()->GetBinLowEdge(nXbins_+1) == xMax_);
  assert(hSrcWeights->GetNbinsY() == nYbins_);
  assert(hSrcWeights->GetYaxis()->GetBinLowEdge(1) == yMin_);
  assert(hSrcWeights->GetYaxis()->GetBinLowEdge(nYbins_+1) == yMax_);

  hSrcWeightsTable_ = hSrcWeights;
  return;
}

TH2F* I3TH2FGenerator::GetSrcWeightsTable() {
  return hSrcWeightsTable_;
}

// This method draws a random (ra,dec) from hSrcWeightsTable_, a TH2F
void I3TH2FGenerator::GenerateCoord() {

  double ranRaDeg=0, ranDecDeg=0;
  hSrcWeightsTable_->GetRandom2(ranRaDeg, ranDecDeg);
  testCoord_ = EquatorialDeg(ranRaDeg, ranDecDeg);

}


int I3TH2FGenerator::GenerateEventEntry() {

  // This method now does everything we need for this location
  SetCandidateEvents();
  
  if (weightMax_<=0.) {
    log_fatal("FATAL: generating source event, but weights are zero!!!\n");
    return -1;  // this is going to crash... what else can you do?
  }

  int entry;
  bool found = false;
  do {
    entry = int ( random_uniform(0.,candidateEvents_.size() ) );
    double prob = random_uniform(0., weightMax_);
    double weight = candidateEvents_[entry].GetMCParams().srcWeight;
    if (prob < weight) { found = true; }
  } while ( !found );  

  return entry;
}


I3Event I3TH2FGenerator::GenerateEvent() {
  int entry = GenerateEventEntry();
  return candidateEvents_[entry];
}
