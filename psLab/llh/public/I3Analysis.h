#ifndef LLH_I3ANALYSIS_H_
#define LLH_I3ANALYSIS_H_

#include <vector>

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/classes.h"
#include "llh/public/EventTimeModule.h"
#include "llh/public/I3Event.h"
#include "llh/public/I3SignalGenerator.h"
#include "llh/public/SimpleEnergyProb.h"


class I3Analysis : public AnalysisSet {
 protected: 

  bool randomizeBase_;
  bool randomizeSrc_;
  bool addSignalToEProb_;
  double baseThinningProb_;

  vector<I3Event> baseEvents_;  // Does Not Change 
  vector<I3Event> modifiableBaseEvents_; // Scrambles; What eList_ points to
  vector<I3Event> sourceEvents_; // also what eList_ points to
  EventPtrList eList_;

  BkgSpaceProb *bkgSpaceProb_;
  SimpleEnergyProb *eProb_; 

  EventTimeModule* evTimeModulePtr_;

  // Protected Functions

  void BindBaseEvents();

  // nSrcEvents<0 ==> generate random number based on source mean
  // nSrcEvents>=0  ==> specifically generate nSrcEvents
  void GenerateSrcEvents(vector<I3Event>& SrcEventVector, int nSrcEvents);

 public:

  I3Analysis(); 
  ~I3Analysis();


  //  *  CONFIGURATION SETTINGS  *  //

  void SetRandomizeBase(bool randomize) { randomizeBase_ = randomize; }
  void SetRandomizeSrc(bool randomize) { randomizeSrc_ = randomize; }
  void SetAddSignalToEProb(bool add) { addSignalToEProb_ = add; }
  void SetBaseThinningProb(double prob) { baseThinningProb_ = prob; }

  void SetBaseEvents(const vector<I3Event> &inputEvents);

  // creates a clone copy which must eventually be deleted
  void SetBkgSpaceProb(const BkgSpaceProb& bsp);

  // Note this will get modified each time the data + injected signal 
  // is passed to eProb_ for updating the bkg tables!
  // (The problem with making a local copy is derived classes..., we 
  // might make a Clone function for EnergyProb's)
  void SetEnergyProb(SimpleEnergyProb &inputEProb);


  // creates a clone copy which must eventually be deleted
  void SetSource(const SourceModule& src);

  // The Module will keep track of event times; note that it is not
  // part of the I3Analysis object itself, so it should *NOT* be deleted
  void SetEventTimeModulePtr( EventTimeModule* evTimeModule);

  // -- This is now handled inside of EventTimeModule
  //  void SetTimeRange(double a, double b);
  
  //  *  OPERATIONS  *  //

  double GetBaseThinningProb() const { return baseThinningProb_; }
  const BkgSpaceProb* GetBkgSpaceProbFn() const { return bkgSpaceProb_; }
  const EnergyProb* GetEnergyProbFn() const { return eProb_; }
  const SourceModule* GetSource() const { return srcModule_; }
  EventTimeModule* GetEventTimeModulePtr() { return evTimeModulePtr_; }

  // Returns current number density, with fake signal added to data
  double BkgNumberDensity(const Coord& coord) const;

  // nSrcEvents<0 ==> generate random number based on source mean
  // nSrcEvents>=0  ==> specifically generate nSrcEvents
  void GenerateDataSet_with_nSrcEvents(int nSrcEvents);

  void GenerateDataSet() { GenerateDataSet_with_nSrcEvents(-1); }

  void UseRealData();

  const EventPtrList* GetEventPtrList() const { return &eList_; }

};


#endif // LLH_I3ANALYSIS_H_

