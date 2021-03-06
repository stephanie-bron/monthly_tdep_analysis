#include "llh/public/I3Analysis.h"

#include "TMath.h"

#include "rootExt/public/log_report.h"
#include "rootExt/public/randomfunctions.h"

#include "llh/public/CoordEquatorialDeg.h"


// PROTECTED FUNCTIONS

void I3Analysis::BindBaseEvents() {
  // wait until all of these ingredients have been set
  if (baseEvents_.size()>0 && bkgSpaceProb_ && eProb_) {
    for (int i=0; i<int(baseEvents_.size()); ++i) {
      baseEvents_[i].SetAnalysisSet(this); 
      baseEvents_[i].SetBkgSpaceProb(bkgSpaceProb_);
      baseEvents_[i].SetEnergyProb(eProb_); //here is the line where the energy pdf is created
    }
    modifiableBaseEvents_ = baseEvents_;
  }   
}


// nSrcEvents<0 ==> generate random number based on source mean
// nSrcEvents>=0  ==> specifically generate nSrcEvents
void I3Analysis::GenerateSrcEvents(vector<I3Event>& SrcEventVector, 
				   int nSrcEvents)
{
  if (nSrcEvents != 0 && !srcModule_) {
    log_fatal("Source was not set... cannot generate source events.\n");
  }

  SrcEventVector.clear();

  int nGeneratedSrcEvents;
  if (nSrcEvents<0) {
    nGeneratedSrcEvents = random_poisson( srcModule_->GetMeanSrcNev() );
  } else {
    nGeneratedSrcEvents = nSrcEvents;
  }
  
  if (nGeneratedSrcEvents){
    I3SignalGenerator* src = dynamic_cast<I3SignalGenerator*>(srcModule_);
    if (!src) {
      log_fatal("ERROR: Expected SourceModule of type I3SignalGenerator.\n");
    }
    
    for (int i=0; i<nGeneratedSrcEvents; ++i) {
      SrcEventVector.push_back( src->GenerateEvent() );
      SrcEventVector[i].SetAnalysisSet(this);
      SrcEventVector[i].SetBkgSpaceProb(bkgSpaceProb_);
      SrcEventVector[i].SetEnergyProb(eProb_);
    }
  }


  /*for (int i=0; i<nGeneratedSrcEvents; ++i) {
    SrcEventVector.push_back( src->GenerateEvent() );
    SrcEventVector[i].SetAnalysisSet(this);
    SrcEventVector[i].SetBkgSpaceProb(bkgSpaceProb_);
    SrcEventVector[i].SetEnergyProb(eProb_);
  }*/
}


// PUBLIC FUNCTIONS


I3Analysis::I3Analysis() :
  randomizeBase_(true), 
  randomizeSrc_(false),
  addSignalToEProb_(true),
  baseThinningProb_(0.),
  bkgSpaceProb_(NULL),
  eProb_(NULL),
  evTimeModulePtr_(NULL)
{ }


I3Analysis::~I3Analysis()
{ 
  if (bkgSpaceProb_) { delete bkgSpaceProb_; }
  if (srcModule_) { delete srcModule_; }
}


void I3Analysis::SetSource(const SourceModule& src) {
  if (srcModule_) { delete srcModule_; }
  // this 'new' copy must be deleted in ~I3Analysis and whenever reset (here) 
  srcModule_ = src.Clone();
}


// Returns current number density, with fake signal added to data
double I3Analysis::BkgNumberDensity(const Coord& coord) const {
  return bkgSpaceProb_->GetBkgProbDensity(coord) * eList_.GetSize();
}


void I3Analysis::SetBaseEvents(const vector<I3Event> &inputEvents)
{
  baseEvents_ = inputEvents;
  BindBaseEvents();
}


void I3Analysis::SetBkgSpaceProb(const BkgSpaceProb& bsp) {
  if (bkgSpaceProb_) { delete bkgSpaceProb_; }
  // this 'new' copy must be deleted in ~I3Analysis and whenever reset (here) 
  bkgSpaceProb_ = bsp.Clone();
  BindBaseEvents();
}


// Note this will get modified each time the data + injected signal 
// is passed to eProb_ for updating the bkg tables!
// (The problem with making a local copy is derived classes..., we 
// might make a Clone function for EnergyProb's)
void I3Analysis::SetEnergyProb(SimpleEnergyProb &inputEProb) { 
  eProb_ = &inputEProb; 
  BindBaseEvents();
}


// The Module will keep track of event times; note that it is not
// part of the I3Analysis object itself, so it should *NOT* be deleted
void I3Analysis::SetEventTimeModulePtr( EventTimeModule* evTimeModule) {
    evTimeModulePtr_ = evTimeModule;
  }


void I3Analysis::GenerateDataSet_with_nSrcEvents(int nSrcEvents) {
  eList_.Clear();
  
  if ( !evTimeModulePtr_ ) {
    evTimeModulePtr_ = new EventTimeModule();

    // (very minor memory leak here if never deleted, but it only happens
    //  once so not very worrisome)

    // times from generic module with default constructor 
    // are only generated over one sidereal day... this is *supposed* to look
    // wierd if someone is expecting a realistic distribution of times
    // during the year...

    // This maintains backward compatibility for time-independent scripts...
    // eventually we may think to remove this and force _all_ scripts to
    // specify how time is handled.
  }

  // Clear the list (if any) of used times which the module tracks internally
  evTimeModulePtr_->ResetUsedTimes();


  // Get Modifiable version of Base (i.e. background data) events

  for (unsigned int i=0; i < modifiableBaseEvents_.size(); ++i) {

    // Use this to sample only a fraction of base 
    if (baseThinningProb_ > 0. && baseThinningProb_ < 1.) {
      if ( random_uniform(0.,1.) >= baseThinningProb_ ) {
	continue;  // skip if random number is ouside sampling range
      }
    }

    // If desired, randomize time and RA, according to method in evTimeModule
    if (randomizeBase_) {
      evTimeModulePtr_->RandomizeEvent( modifiableBaseEvents_[i] );
    }

    eList_.AddEvent(&(modifiableBaseEvents_[i]));
  }
  

  // Get Source events
  
  sourceEvents_.clear();
  EventPtrList sList;

  if(nSrcEvents != 0)
    {
      GenerateSrcEvents(sourceEvents_, nSrcEvents);

      for (unsigned int i=0; i < sourceEvents_.size(); ++i) {
	// If desired, scramble RA for this event
	if (randomizeSrc_) {
	  evTimeModulePtr_->RandomizeEvent( sourceEvents_[i] );
	}
	
	sList.AddEvent( &(sourceEvents_[i]));
	eList_.AddEvent(&(sourceEvents_[i]));
      }
    }
  // are there cases where this is really not necessary, and could speed things
  // up if the sorting were skipped?  (and would it matter much?)
  eList_.SortByTime();


  if ( sList.GetSize() > 0 ) {
    bkgSpaceProb_->FixToBasePlusEvents(sList);
  } 
  else {
    bkgSpaceProb_->FixToBase();
  }

  if (addSignalToEProb_) {
    eProb_->SetTableBkg(eList_);
  }
  
}


void I3Analysis::UseRealData() {
  bkgSpaceProb_->FixToBase();
  modifiableBaseEvents_ = baseEvents_;
  
  eList_.Clear();
  for (unsigned int i=0; i<modifiableBaseEvents_.size(); ++i) {
    eList_.AddEvent( &(modifiableBaseEvents_[i]) );
  }

  // are there cases where this is really not necessary, and could speed things
  // up if the sorting were skipped?  (and would it matter much?)
  eList_.SortByTime();

  eProb_->SetTableBkg(eList_);
}
