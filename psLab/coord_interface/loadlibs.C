
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "coord_interface";

  vector<string> project_dependencies;
  // none

  // For now, only thing this project does is load coordinate-service
  // icetray library and provide public headers for interface

  BasicLoad(this_project, project_dependencies);
}
