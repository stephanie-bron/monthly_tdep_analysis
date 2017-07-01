
{

  TString cutName = "CutA2";

  //TCut testCut = "(s16g1blogl+s16g2blogl)-s32logl>35 && (s16t1blogl+s16t2blogl)-s32logl>35 && mLDirC>=200 && mNDirC>=5";


  gROOT->ProcessLine(".x loadlibs.C");

  gROOT->ProcessLine(".x macro_loadCutA2.C");

  bool SUPER_MACRO_USE_ENERGY;

  bool SUPER_MACRO_DISCOVERY = true;

  TString outdir = "TestCuts/";
  gSystem->mkdir(outdir);
  TString filebase = "ic40_zen_";
  TString filebaseSens = "Neyman_sin_Mue_MPE_SigmaCor_";
  TString filebaseDisc = "Disc_sin_Mue_MPE_SigmaCor_";
  TString filename;


  // Sensitivities //
  spectralIndex = -2;
  filename = filebase + "Em2_" + filebaseSens + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_sensitivity_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");

  spectralIndex = -1.5;
  filename = filebase + "Em1.5_" + filebaseSens + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_sensitivity_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");

  // Discovery Potentials //
  spectralIndex = -2;
  filename = filebase + "Em2_" + filebaseDisc + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");

  spectralIndex = -1.5;
  filename = filebase + "Em1.5_" + filebaseDisc + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");

  spectralIndex = -3;
  filename = filebase + "Em3_" + filebaseDisc + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");

  // Out of Order
  spectralIndex = -3;
  filename = filebase + "Em3_" + filebaseSens + cutName + ".root";
  cout << "\nStarting for Filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_sensitivity_all_zenith.C");
  cout << "Writing filename: " << filename << endl;
  gROOT->ProcessLine(".x macro_discovery_zenith_save.C(\""+outdir+filename+"\")");


}
