TString LOADTREE_86_II  = "/data/user/tcarver/Data/PS/root_files/IC86_2012/";
TString LOADTREE_86_III = "/data/user/tcarver/Data/PS/root_files/IC86_2013/";
TString LOADTREE_86_IV  = "/data/user/tcarver/Data/PS/root_files/IC86_2014/";

TTree* LoadTree_IC86_II_III_IV_nugen_numu_E_1() {
    vector <TString> nufiles;
    nufiles.push_back("Merged_11069_downgoing_6089filesIC86.2012.root");
    nufiles.push_back("Merged_11069_upgoing_6089filesIC86.2012.root");
    nufiles.push_back("Merged_11070_downgoing_4829filesIC86.2012.root");
    nufiles.push_back("Merged_11070_upgoing_4829filesIC86.2012.root");
    return LoadTree_IC86_II_III_IV_nugen_numu(nufiles);
}
  
TTree* LoadTree_IC86_II_III_IV_nugen_numu(vector <TString> nufiles) {
    TChain *tree  =new TChain("MasterTree");

    for (unsigned int j=0;j<nufiles.size();j++){
        tree->Add(LOADTREE_86_II+nufiles[j]);
        std::cout << "LOADTREE_86_II+nufiles[j]" << LOADTREE_86_II+nufiles[j] << std::endl;
    }
    
    tree->GetEntries();
    tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon.energy");
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
    tree->SetAlias("NChan",                 "HitMultiplicityValuesIC.n_hit_doms");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
    tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
    tree->SetAlias("mcPrimary_Zenith_rad",  "MCPrimary1.zenith");
    tree->SetAlias("mcPrimary_Azimuth_rad", "MCPrimary1.azimuth");
    tree->SetAlias("mcPrimary_Energy_GeV",  "MCPrimary1.energy");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",               "I3EventHeader.Event");
    tree->SetAlias("mcOneWeight", "I3MCWeightDict.OneWeight / (I3MCWeightDict.NEvents*((TMath::FloorNint(RunID/100000)==11069)*6089 + (TMath::FloorNint(RunID/100000)==11070)*4829))");

    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"DelAng", "acos( sin(MCPrimary1.zenith)*sin("+TString(recos[i])+".zenith)*cos(MCPrimary1.azimuth-"+TString(recos[i])+".azimuth) + cos(MCPrimary1.zenith)*cos("+TString(recos[i])+".zenith))*TMath::RadToDeg()");
    }
  
    tree->SetAlias("mcTotalGeneratedEvents", "1.*1.");
    return tree;
}

TTree* LoadTree_IC86_II_III_IV_GoodRuns_Full() {
    TString upfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_upgoing.root";
    std::cout << "upfiles2012 " << upfiles2012 << std::endl;
    TString downfiles2012 = LOADTREE_86_II+"Merged_IC86.2012_downgoing.root";
    std::cout << "downfiles2012 " << downfiles2012 << std::endl;
    TString upfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_upgoing.root";
    std::cout << "upfiles2013 " << upfiles2013 << std::endl;
    TString downfiles2013 = LOADTREE_86_III+"Merged_IC86.2013_downgoing.root";
    std::cout << "downfiles2013 " << downfiles2013 << std::endl;
    TString upfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_upgoing.root";
    std::cout << "upfiles2014 " << upfiles2014 << std::endl;
    TString downfiles2014 = LOADTREE_86_IV+"Merged_IC86.2014_downgoing.root";  
    std::cout << "downfiles2014 " << downfiles2014 << std::endl;

    TChain *tree = new TChain("MasterTree");
    tree->Add(downfiles2012);
    tree->Add(upfiles2012);
    tree->Add(downfiles2013);
    tree->Add(upfiles2013);
    tree->Add(downfiles2014);
    tree->Add(upfiles2014);

    tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

    //here are the infos about on the dates of the runs
    TString livetimeTotalStr = " 1. * 91371456.0"; // livetime in sec (II+III = 59644512.0)
//    TString tmin = " 1. * 56062.4207"; //in MJD ->  '2012-05-15 10:05:48.480' in normal date
//    TString tmax = " 1. * 57160.0410"; // for III = 56783.5781"; 57160.0410 = '2015-05-18 00:59:02.400'
    // set wrong times just in order to go until the date we want for the sources of interest
    TString tmin = " 1. * 56592.4"; //in MJD ->  '2013-10-27 09:36:00.000' in normal date
    TString tmax = " 1. * 57314.0"; // '2015-10-19 00:00:00.000'

    // SET META-INFORMATION:
    // These aliases store meta-information for the tree
    tree->SetAlias("livetimeTotal", livetimeTotalStr);
    tree->SetAlias("tmin", tmin);
    tree->SetAlias("tmax", tmax);
  
    //--------------------//
    tree->SetAlias("timeMJD","timeMJD.value");
    // OTHER, EXTRA ALIASES:
  
    tree->SetAlias("truncEn",          "SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon.energy");
    tree->SetAlias("muexEn",                "SplineMPEMuEXDifferential.energy");
    tree->SetAlias("NChan",                 "HitMultiplicityValuesIC.n_hit_doms");
    string recos[2]={"SplineMPE","MuEXAngular4"};
    for (int i=0; i<2; i++) {
        tree->SetAlias(TString(recos[i])+"Zd",              TString(recos[i])+".zenith*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Ad",              TString(recos[i])+".azimuth*TMath::RadToDeg()");
        tree->SetAlias(TString(recos[i])+"Zr",              TString(recos[i])+".zenith");
        tree->SetAlias(TString(recos[i])+"Ar",              TString(recos[i])+".azimuth");
    }
    tree->SetAlias("mpbSigmaDeg",           "ParabSigma.value*TMath::RadToDeg()");
    tree->SetAlias("mMuEXSigmaDeg",         "MuEXAngular4_Sigma.value*TMath::RadToDeg()");
    tree->SetAlias("CRSigmaDeg",            "CRSigmaDeg.value");
    tree->SetAlias("SplinePbSigmaDeg",      "sqrt(pow(SplineMPEParaboloidFitParams.err1,2)+pow(SplineMPEParaboloidFitParams.err2,2))/sqrt(2)*TMath::RadToDeg()");
    tree->SetAlias("RunID",                 "I3EventHeader.Run");
    tree->SetAlias("EventID",                 "I3EventHeader.Event");
    return tree;
}
