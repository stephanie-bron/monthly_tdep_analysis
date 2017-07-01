
{

  // THIS PROJECT

  TString projectSrcDir = "llhTimeDep/public";

  // LIST OF LIBRARIES

  vector<string> buildList;

  buildList.push_back("TimePdfCollection.h");
  buildList.push_back("BlockLevel.h");
  buildList.push_back("FracUpTime.C");
  buildList.push_back("LocalCoordBkgProb.C");
  buildList.push_back("NewLlhBlockTime.C");
  buildList.push_back("NewLlhGausTime.C");
  buildList.push_back("NewLlhPeriodicTime.C");
  buildList.push_back("GetSrcDeltaEterm.C");
  buildList.push_back("MultiBlockAnalysisFn.C");
  buildList.push_back("NewLlhBoxTime.C");
  //buildList.push_back("MultiTDAnalysisFn.C");
  buildList.push_back("MultiPeriodicAnalysisFn.C");
  
  // For using vectors in CINT macros:
  // (Put any new classes which you need from this project in here
  //  buildList.push_back("VectorLoader.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
