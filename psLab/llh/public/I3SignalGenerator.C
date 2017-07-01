#include "llh/public/I3SignalGenerator.h"

#include "rootExt/public/randomfunctions.h"
#include "fluxus/public/FluxFunction.h"


//
// I3MultiSignalGenerator
//



// Need a deep copy, since new copies must be made for signalPtr's
I3SignalGenerator* I3MultiSignalGenerator::Clone() const {
  I3MultiSignalGenerator* newMultiPtr = new I3MultiSignalGenerator();
  newMultiPtr->livetime_ = livetime_;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    newMultiPtr->AddSignal(*(signalPtrVect_[i]), enhanceFactorVect_[i]);
  }
  return newMultiPtr;
}


void I3MultiSignalGenerator::SetLivetime(double livetime) { // THIS IS FRAUGHT
  livetime_ = livetime;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    signalPtrVect_[i]->SetLivetime(livetime_);
  }
}

void I3MultiSignalGenerator::SetTimePdf(TimePdf * tPdf) {
  timePdf_ = tPdf;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    signalPtrVect_[i]->SetTimePdf(tPdf);
  }
  
  double factor;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {  
    factor = signalPtrVect_[i]->GetTimePdf()->GetNorm();
    SetEnhanceFactor(i, factor);
  }
  
  // OKAY, above is a first trial at adding time to the MultiSignalGenerator
  // ... Chad added in this magical Enhancement Factor, and the fraction
  // of the timePdf in a given dataset seems like a good use of that.
  // I think that will break / obfuscate other things, like the conversion
  // from nevs to a comparable steady flux (or fluence, whenever I get to 
  // that). So I'm making some "Base" functions which neglect the Enhance
  // factor as a quick way to check with previous (pre-multi) things.
  
  // Or should it be the opposite and the regular should be "Time"? meh.
  
  // Another problem: the data knows about tmin and tmax, not the source
  // for now (that used to be the timePdf's job which was set by hand).
  // This doesn't do the CheckTimeBounds(t1,t2) to properly set the
  // normalization yet.
  
  // This is set up in the MultiArk
  
}

void I3MultiSignalGenerator::AddSignal(const I3SignalGenerator& signal, 
					  double enhanceFactor) 
{
  // This makes a private copy of the added signal; must delete in destructor
  I3SignalGenerator* newSignalPtr = signal.Clone();
  signalPtrVect_.push_back(newSignalPtr);
  enhanceFactorVect_.push_back(enhanceFactor);
}

double I3MultiSignalGenerator::GetMeanSrcNevBase() const { // HMMM
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += signalPtrVect_[i]->GetMeanSrcNev();
  }
  return nev;
}

double I3MultiSignalGenerator::GetMeanSrcNev() const {
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += enhanceFactorVect_[i] * signalPtrVect_[i]->GetMeanSrcNev();
  }
  return nev;
}




double I3MultiSignalGenerator::
GetMeanSrcNevForFluxModelBase(const FluxBase& fluxModel) const { // HMMM
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += signalPtrVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return nev;
}

double I3MultiSignalGenerator::
GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += enhanceFactorVect_[i] * 
      signalPtrVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return nev;
}





I3Event I3MultiSignalGenerator::GenerateEvent() {
  double meanNevSum = GetMeanSrcNev();
  double ranNum = random_uniform(0., meanNevSum);

  // start adding up the nev (as in GetMeanSrcNev), and stop when 
  // random number is less than sum

  double sum = 0.;
  unsigned int i;
  for (i = 0; i< signalPtrVect_.size(); ++i) {
    sum += enhanceFactorVect_[i] * signalPtrVect_[i]->GetMeanSrcNev();
    if (ranNum <= sum) {
      break;
    }
  }

  // 'i' is now the signal which we want to use:
  return signalPtrVect_[i]->GenerateEvent();
}




//
// I3PointGenerator
//



I3PointGenerator::I3PointGenerator() : 
  livetime_(0),
  weightSum_(0),
  weightMax_(0),
  timeAzBins_(1),
  sourceTimePdf_(NULL)  
{ }


I3PointGenerator::I3PointGenerator(
				   const vector<I3Event>& inputEvents, 
				   const FluxBase& fluxModel,
				   const EquatorialDeg& sourceCoord, 
				   double livetime) :
  livetime_(livetime)
{
  SetCandidateEvents(inputEvents, fluxModel, sourceCoord);
  SetTimeAzBins(1);
  sourceTimePdf_ = NULL;
}

I3PointGenerator::I3PointGenerator(
				   const vector<I3Event>& inputEvents, 
				   const FluxBase& fluxModel,
				   const EquatorialDeg& sourceCoord, 
				   double livetime,
 				   TimePdf * tPdf) :
  livetime_(livetime)
{
  SetCandidateEvents(inputEvents, fluxModel, sourceCoord);
  SetTimeAzBins(1);
  SetTimePdf(tPdf);
}

void I3PointGenerator::SetCandidateEvents(const vector<I3Event>& inputEvents,
					  const FluxBase& fluxModel,
					  const EquatorialDeg& sourceCoord) { 
  sourceCoord_ = sourceCoord;
  candidateEvents_.clear();
  candidateEvents_ = inputEvents;
  SetCandidateWeights(fluxModel);
}


// Weights do not include the livetime (therefore, livetime can be changed 
// independently at any time).  Effectively, this means each weight stands
// for "the mean number of such events in 1 sec"
// (since sec is unit of livetime).

void I3PointGenerator::SetCandidateWeights(const FluxBase& fluxModel) {
  weightSum_ = 0.;
  weightMax_ = 0.;

  for (vector<I3Event>::iterator e = candidateEvents_.begin();
       e != candidateEvents_.end();
       e++) 
  {
    I3MCParameters mc = e->GetMCParams();

    double weight = fluxModel.GetFlux(mc.mcEnergy) * mc.PS_FlatSpectrumRate;
    // this gets multiplied by livetime for number of events, so
    // effectively it could be considered "number of events for 1 sec" here.

    if (weight<0) { 
      log_error("ERROR: Tried to assign weight<0\n"); 
      weight = 0.;
    }
      
    mc.srcWeight = weight;
    e->SetMCParams(mc);

    weightSum_ += weight;
    if (weight > weightMax_) 
      {weightMax_ = weight;}
  }
}


double I3PointGenerator::
GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double sum = 0.;

  for (int i=0; i<int(candidateEvents_.size()); ++i) {
    const I3Event& e = candidateEvents_[i];
    I3MCParameters mc = e.GetMCParams();
    double weight = fluxModel.GetFlux(mc.mcEnergy) * mc.PS_FlatSpectrumRate;
    // this gets multiplied by livetime for number of events, so
    // effectively it could be considered "number of events for 1 sec" here.

    if (weight<0) { 
      log_error("ERROR: Found negative weight!\n"); 
      weight = 0.;
    }
    sum += weight;
  }
  return sum * livetime_;
}

/*
I3Event I3PointGenerator::GenerateEvent() {
  I3Event *e = NULL;
  Time t = sourceTimePdf_->GenerateEventTime();
  double srcAzDeg =  145.6453 + fmod( t.GetMJD(), 0.997269566)/0.997269566*360. - sourceCoord_.GetRa();
  while (srcAzDeg>360) srcAzDeg -= 360.;
  if (srcAzDeg<0) srcAzDeg += 360.;

  int azBin, azBinE; // this checks that the event is in the right part of the local coordinate space
                     // to be consistent with arriving at the signal event time.
  if (timeAzBins_) {
    azBin = (int) srcAzDeg/(360/timeAzBins_);
  } else {
    azBin=0;
  }

  do {
    int randomID = int ( random_uniform(0.,candidateEvents_.size() ) );
    double prob = random_uniform(0., weightMax_);
    double weight = candidateEvents_[randomID].GetMCParams().srcWeight;
    if (timeAzBins_) {
      azBinE = (int) candidateEvents_[randomID].GetParams().recoAzimuthDeg/(360/timeAzBins_);
    } else {
      azBinE = 0;
    }

    if ( prob < weight && azBin==azBinE ) {
      e = &candidateEvents_[randomID];
    }
  } while (!e);

  e->SetTime( t );
  return *e;

}*/


int I3PointGenerator::GenerateEventEntry() {
  if (weightMax_<=0.) {
    log_fatal("FATAL: generating source event, but weights are zero!!!\n");
    return -1;  // this is going to crash... what else can you do?
  }
  
  double srcAzDeg=0.;
  Time t = Time(0.);
  if ( sourceTimePdf_ ) {
    t = sourceTimePdf_->GenerateEventTime();
  
    srcAzDeg =  145.6453 + fmod( t.GetMJD(), 0.997269566)/0.997269566*360. - 
                       sourceCoord_.GetRa();
                           
    while (srcAzDeg>360) srcAzDeg -= 360.;
    if (srcAzDeg<0) srcAzDeg += 360.;
    
  } 

  int azBin, azBinE; // this checks that the event is in the right 
                     // part of the local coordinate space
                     // to be consistent with arriving at the signal event time.
  if (timeAzBins_) {
    azBin = (int) srcAzDeg/(360./timeAzBins_);
  } else {
    azBin=0;
  }

  int entry;
  bool found = false;
  do {
    entry = int ( random_uniform(0.,candidateEvents_.size() ) );
    double prob = random_uniform(0., weightMax_);
    double weight = candidateEvents_[entry].GetMCParams().srcWeight;
    
    if (timeAzBins_) {
      azBinE = (int) candidateEvents_[entry].GetParams().recoAzimuthDeg/(360./timeAzBins_);
    } else {
      azBinE = 0;
    }
        
    if (prob < weight && azBin==azBinE ) { found = true; }
  } while ( !found );  

  candidateEvents_[entry].SetTime( t );

  return entry;
}


I3Event I3PointGenerator::GenerateEvent() {
  int entry = GenerateEventEntry(); // Time has to be set when we get the entry number
  return candidateEvents_[entry];
}


