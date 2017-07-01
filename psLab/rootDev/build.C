
{

  // THIS PROJECT

  TString projectSrcDir = "rootDev/public";

  // LIST OF LIBRARIES

  vector<string> buildList;

  buildList.push_back("SetPad.h");
  //  buildList.push_back("ColorChooser.C");
  buildList.push_back("CreatePalette.C");
  buildList.push_back("EventController.C");
  buildList.push_back("tempMH.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
