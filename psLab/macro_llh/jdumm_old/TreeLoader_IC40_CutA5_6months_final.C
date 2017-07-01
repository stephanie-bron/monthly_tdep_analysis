
TTree* LoadTree_febIC40_CutA5_nugen1882_small() {

  char *files = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA5/CutA5_Nugen1882_01000_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(files);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN


  int nFiles = 1000; // Only 869 files right now
  int generatedEventsPerFile = 5000; // NuMu + NuMuBar per file
  int mcTotalGeneratedEvents = nFiles * generatedEventsPerFile;


  // SET HANDLES:
  // These aliases provide consistent handles across different 
  // types root files (analysis-tree, flat-ntuple, simple-extractor),
  // so set accordingly:

  tree->SetAlias("mcPrimary_Zenith_rad",  "mcZr");
  tree->SetAlias("mcPrimary_Azimuth_rad", "mcAr");
  tree->SetAlias("mcPrimary_Energy_GeV",  "mcEn");
  tree->SetAlias("mcOneWeight",           "OW");
  tree->SetAlias("mueEn",                 "mmueEn");
  tree->SetAlias("NChan",                 "NCh");

  tree->SetAlias("z",                     "cos(mZd*TMath::DegToRad())");

  tree->SetAlias("BartolFluxWeight",      "BartolFluxWeightForOneFile");

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("mcTotalGeneratedEvents", TStringify(mcTotalGeneratedEvents));


  //--------------------//


  // OTHER, EXTRA ALIASES:

  tree->SetAlias("s32Zr","s32Zd*TMath::DegToRad()");
  tree->SetAlias("s32Ar","s32Ad*TMath::DegToRad()");
  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");


  return tree;
}



TTree* LoadTree_febIC40_CutA5_data6months_GoodRuns_small() {

   char *filespec10 = "/net/user/jdumm/std-processing/releases/V03-02-00/release64/CutA5/CutA5_6months_GoodRuns_small.root";
  TChain *tree = new TChain("Level2");
  tree->Add(filespec10);
  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  double livetimeTotal = 175.55 * 86400.; // seconds

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree

  tree->SetAlias("livetimeTotal", TStringify(livetimeTotal));

  //--------------------//

  // OTHER, EXTRA ALIASES:
  tree->SetAlias("timeMJD","MJDay+(MJSec+MJNanoSec*1e-9)/86400.");

  tree->SetAlias("NChan","NCh");
  tree->SetAlias("mueEn","mmueEn");
  tree->SetAlias("s32Zr","s32Zd*TMath::DegToRad()");
  tree->SetAlias("s32Ar","s32Ad*TMath::DegToRad()");
  tree->SetAlias("mZr","mZd*TMath::DegToRad()");
  tree->SetAlias("mAr","mAd*TMath::DegToRad()");

  tree->SetAlias("z","cos(mZd*TMath::DegToRad())");

  return tree;
}
