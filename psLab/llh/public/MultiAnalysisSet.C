#include "llh/public/MultiAnalysisSet.h"

#include "rootExt/public/randomfunctions.h"

#include "llh/public/TimePdf.h"
#include "llh/public/I3Event.h"
#include "llh/public/I3SignalGenerator.h"
#include "llh/public/I3Analysis.h"




void MultiAnalysisSet::ConstructEventPtrList() {
  evPtrListMulti_.Clear();
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    const EventPtrList* tempList = aSetVect_[i]->GetEventPtrList();
    for (int j=0; j<tempList->GetSize(); ++j) {
      evPtrListMulti_.AddEvent(tempList->GetEvent(j));
    }
  }
}


void MultiAnalysisSet::AddAnalysisSet(AnalysisSet* aSet) { 
  aSetVect_.push_back(aSet); 
  ConstructEventPtrList();
}


double MultiAnalysisSet::BkgNumberDensity(const Coord& coord) const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->BkgNumberDensity(coord);
  }
  return sum;
}


double MultiAnalysisSet::GetMeanSrcNev() const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->GetMeanSrcNev();
  }
  return sum;
}

double MultiAnalysisSet::GetMeanSrcNevTime() const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->GetSource()->GetTimePdf()->GetNorm() * aSetVect_[i]->GetMeanSrcNev();
  }        // Woof!
  return sum;
}

double MultiAnalysisSet::
GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return sum;
}


double MultiAnalysisSet::GetEnergyQuantile(double prob) const 
{
  int nSets = aSetVect_.size();
  
  int nbins = 100;

  double lowedge = 2.;
  double highedge = 10.;
  
  TH1D * srcHisto = new TH1D("srcHisto","srcHisto", nbins, lowedge, highedge);
  
  for (int i = 0; i < nSets; ++i)
    {
      
      const SourceModule *srcModule = aSetVect_[i]->GetSource();
      
      const I3PointGenerator *i3point = dynamic_cast<const I3PointGenerator*>(srcModule);
	
      
      for (int j=0; j<int(i3point->candidateEvents_.size()); ++j) {
	const I3Event& e = i3point->candidateEvents_[j];
	I3MCParameters mc = e.GetMCParams();
	
	
	if(log10(mc.mcEnergy) < lowedge || log10(mc.mcEnergy) > highedge)
	  {
	    log_warn("Energy out histogram range: %f\n", log10(mc.mcEnergy));
	  }
	
	srcHisto->Fill(log10(mc.mcEnergy),mc.srcWeight);   
      }
      
    }
  
  double quantile;
  
  srcHisto->GetQuantiles(1, &quantile, &prob);
  
  delete srcHisto;
  
  return pow(10., quantile);
    
}


double MultiAnalysisSet::GetFluxScaleForNev(double nev) const
{
  // Use contribution from a single, non-zero set to determine flux
  // (This assumes that same flux was used for all sets... 
  //  otherwise this question is nonsensical to begin with)

  double totalMeanSrcNev = GetMeanSrcNev();
  int nSets = aSetVect_.size();
  for (int i=0; i<nSets; ++i) {
    double partialMeanSrcNev = aSetVect_[i]->GetMeanSrcNev();

    // check that this contribution is non-zero and non-negligible
    if (partialMeanSrcNev > totalMeanSrcNev/(nSets*2.)) {
      // at least one set has to pass this test, unless all contribute zero

      double thisSetNev = nev * (partialMeanSrcNev / totalMeanSrcNev);
      return aSetVect_[i]->GetFluxScaleForNev(thisSetNev);
    }
  }
  return 0.;  // evidently, sources are not generating any signal events
}


void MultiAnalysisSet::GenerateDataSet() {
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    aSetVect_[i]->GenerateDataSet();
  }
  ConstructEventPtrList();
}


void MultiAnalysisSet::GenerateDataSet_with_nSrcEvents(int nSrcEvents) {

  const int nSets = aSetVect_.size();
  vector<int> nSrcPerSet(nSets,0); // initialize all values to zero

  double meanSrcNev_Sum = GetMeanSrcNev();
  double meanSrcNev_TimeSum = 0.;
  
  if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
    meanSrcNev_TimeSum = GetMeanSrcNevTime();  
  }
   
  // randomly decide which DataSet to generate each signal event in,
  // according to the relative weight of meanSrcEvents from each set  
  if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
    // This is the parallel way to get the nsrc if there is a TimePdf set (null by default)
    for (int nev=0; nev<nSrcEvents; ++nev) {
      double sumT = 0.;
      double ranNumT = random_uniform(0., meanSrcNev_TimeSum);
      for (int set=0; set<nSets; ++set) {
        sumT += aSetVect_[set]->GetSource()->GetTimePdf()->GetNorm() * aSetVect_[set]->GetMeanSrcNev();
        if (ranNumT <= sumT) {
  	  // this is the dataset which will get this src event
	      nSrcPerSet[set] += 1;
          break;
        }
      }
    }  
  } else {
    for (int nev=0; nev<nSrcEvents; ++nev) {
      double sum = 0.;
      double ranNum = random_uniform(0., meanSrcNev_Sum);
      for (int set=0; set<nSets; ++set) {
        sum += aSetVect_[set]->GetMeanSrcNev();
        if (ranNum <= sum) {
	  // this is the dataset which will get this src event
      nSrcPerSet[set] += 1;
	  break;
        }
      }
    }
  }
  
  //for (int set=0; set<nSets; ++set) { cout << nSrcPerSet[set] << " " << flush; }
  //cout << endl;

  for (int set=0; set<nSets; ++set) {
    aSetVect_[set]->GenerateDataSet_with_nSrcEvents(nSrcPerSet[set]);
  }
  
  ConstructEventPtrList();
  
}


void MultiAnalysisSet::UseRealData() {
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    aSetVect_[i]->UseRealData();
  }
  ConstructEventPtrList();
}
