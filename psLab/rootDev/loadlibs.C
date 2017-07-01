
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "rootDev";

  vector<string> project_dependencies;
  project_dependencies.push_back("rootExt");

  BasicLoad(this_project, project_dependencies);
}
