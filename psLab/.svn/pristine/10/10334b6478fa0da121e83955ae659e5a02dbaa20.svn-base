#ifndef LLH_MULTIANALYSISSET_H_
#define LLH_MULTIANALYSISSET_H_

#include <vector>

#include "llh/public/classes.h"


class MultiAnalysisSet : public AnalysisSet {
 protected:
  vector<AnalysisSet*> aSetVect_;
  EventPtrList evPtrListMulti_;
  
  // Execute this after actions which change event lists
  void ConstructEventPtrList();

 public:

  MultiAnalysisSet() { }
  ~MultiAnalysisSet() { }

  //  *  CONFIGURATION SETTINGS  *  //

  void AddAnalysisSet(AnalysisSet* aSet);

  //  *  OPERATIONS  *  //

  vector<AnalysisSet*>& GetAnalysisSetVect() { return aSetVect_; }

  double BkgNumberDensity(const Coord& coord) const;

  //There is not one SourceModule defined for the MultiAnalysisSet
  const SourceModule* GetSource() const { return NULL; }

  // THIS WILL HAVE TO BE IMPLEMENTED FOR TIME-DEP CODE!
  EventTimeModule* GetEventTimeModulePtr() { return NULL; }

  double GetMeanSrcNev() const;
  double GetMeanSrcNevTime() const;

  
  double GetEnergyQuantile(double prob) const;
  
  double GetSrcEmin(){ return GetEnergyQuantile(0.05);}
  double GetSrcEmax(){ return GetEnergyQuantile(0.95);}



  double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;


  double GetFluxScaleForNev(double) const;

  void GenerateDataSet();

  void GenerateDataSet_with_nSrcEvents(int nSrcEvents);

  void UseRealData();

  const EventPtrList* GetEventPtrList() const { return &evPtrListMulti_; }
};


#endif // LLH_MULTIANALYSISSET_H_
