
{
  /*
  if (1) {
    int ranSeed = -55;  // set the random seed using a negative integer
    cout << "Using random seed " << ranSeed << "\n\n";
    initialize_ran1(ranSeed);
  }

  random_uniform(0,1); // if ranSeed wasn't set, clock time is used
  */


  // For now, this is where the meta-TTree loading scripts are:
  //gROOT->ProcessLine(".L LoadTreeFns_IC22.C");
  gROOT->ProcessLine(".L TreeLoader_IC40_CutA5_6months_final.C");


  // Globals: filled in this script, and later referred to by other scripts
  double livetime = 0;

  TTree *baseTree;
  TTree *srcTree;

  EventLoader evLoader;



  //
  // CONFIGURE EVENTLOADER WITH BACKGROUND / REAL DATA SAMPLE
  //

  cout << "Configuring Base Data Sample...\n";

  /* Example 1: Use simulation for bkg data */
  if (0) {
    baseTree = LoadTree_mayIC22_nugen651_extuple();
    livetime = 275.7 * 86400.;

    double rateAdjustment = 5114. / 4642.; //  real data / atmNu sim. only
    TString weightString = 
      "BartolFluxWeight*" + TStringify(livetime * rateAdjustment);
    evLoader.SetBkgLoadMethod_PoissonSamplePlus(weightString);
    cout << "\n Simulation rate adjustment: " << rateAdjustment << "\n";

    evLoader.SetTimeMethod_Scramble();
    evLoader.SetZenithRangeDeg(85.,180.,10.);
    // Parameters: zenMin, zenMax, zenithRanShift:
    //   when an event is loaded more than once, its zenith will be randomized
    //   +/- zenRanShift (following phase space on sphere),
    //   but randomization will be confined within Min and Max values given
  }


  /*  Example 2: Use a short, real data sample for bkg data */
  if (0) {
    baseTree = LoadTree_mayIC22_data30day_extuple();
    double dataFileLivetime = GetValueFromTree(baseTree,"livetimeTotal");
    livetime = 275.70 * 86400.;  // THIS IS THE TIME WE WANT TO SCALE UP TO

    double upScaleFactor = livetime / dataFileLivetime;
    evLoader.SetBkgLoadMethod_PoissonSamplePlus(TStringify(upScaleFactor));
    cout << "\n\n Factor " << upScaleFactor << " UPSCALING!\n";

    evLoader.SetTimeMethod_Scramble();
    evLoader.SetZenithRangeDeg(85., 180., 10.);
    // Parameters: zenMin, zenMax, zenithRanShift:
    //   when an event is loaded more than once, its zenith will be randomized
    //   +/- zenRanShift (following phase space on sphere),
    //   but randomization will be confined within Min and Max values given
  }


  /* Example 3: Use the full, real data sample */
  if (1) {
    baseTree = LoadTree_febIC40_CutA5_data6months_GoodRuns_small();
    livetime = GetValueFromTree(baseTree,"livetimeTotal");

    evLoader.SetBkgLoadMethod_Exact();
    
    evLoader.SetTimeMethod_Scramble();  // for blindness!
    //evLoader.SetTimeMethod_Actual("timeMJD"); 
    //    cout << "\n\nUSING ACTUAL TIMES!!!\n\n";
  }

  /*  Use a short, real data sample for bkg data */
  if (0) {
    baseTree = LoadTree_febIC40_CutA5_data161day_small();
    double dataFileLivetime = GetValueFromTree(baseTree,"livetimeTotal");
    livetime = 330. * 86400.;  // THIS IS THE TIME WE WANT TO SCALE UP TO

    double upScaleFactor = livetime / dataFileLivetime;
 cout << "About to Set bkg load method\n";
 cout << TStringify(upScaleFactor) << endl;
    evLoader.SetBkgLoadMethod_PoissonSamplePlus(TStringify(upScaleFactor));
    cout << "\n\n Factor " << upScaleFactor << " UPSCALING!\n";

    evLoader.SetTimeMethod_Scramble();
    //evLoader.SetZenithRangeDeg(0., 40., 10.);
    evLoader.SetZenithRangeDeg(0., 180., 2.0); // Not much wiggle near horizon!
    // Zenith range containing psuedo-events with random-shifted zeniths
  }



  cout << "Livetime (Days): " << livetime/86400. << "\n\n";
  evLoader.SetBkgTree(baseTree);



  //
  // CONFIGURE EVENTLOADER WITH SIGNAL SIMULATION DATA SAMPLE
  //

  // evLoader will load signal events for any specified declination

  double sourceZenWidthDeg;
  // This is the +/- zenith range in degrees to select MC events from,
  // and then rotate to desired source location.
  // (Smaller range is more accurate, but trade-off is lower statistics
  // for signal simulation.  You have to pick something appropriate for the 
  // statistics you have.)

  cout << "Configuring Source Event Sample...\n";

  // Example 1:  E^-1 nuGen file
  if (0) {
    //    srcTree = LoadTree_mayIC22_nugen651_extuple();
    srcTree = LoadTree_finalCuts_mayIC22_nugen651_extuple();
    sourceZenWidthDeg = 4.;
  }

  // Example 1:  E^-1 nuGen file
  if (1) {
    //    DataSet *dsSrc = Load_mayIC22_nugen651_extuple();
    // TTree *srcTree = dsSrc->tree;
    //srcTree = LoadTree_febIC40_L2_nugen1794_small();
    srcTree = LoadTree_febIC40_CutA5_nugen1882_small();
    sourceZenWidthDeg = 0.5;
  }

  evLoader.SetSourceTree(srcTree);
  evLoader.SetSourceZenWidthDeg(sourceZenWidthDeg);



  //
  // SPECIFY CUTS
  //

  if (0) {
    // IC22 PS Cut
    TCut testCut = "g32Zd>80 && pf32PbfStatus==0 && pf32SigmaDeg<3 && !(g32Rlogl>7.8 && g32NdirC<7) && !(g32Rlogl>8.5 && g32NdirC<8) && g32umLogl-g32Logl>15 && sZenMinDeg>70 && b32Logl-g32Logl>30 && g32Rlogl<9.5";
  }
  if (1) {
    TCut testCut = "1";
    cout << "\n   *** NOT APPLYING CUTS ...";
    cout << "   ASSUMING THEY WERE ARE APPLIED TO FILES !!!\n";
  }

  // Cuts can be added or reset any time, and events re-loaded with new cuts
  // (cuts are also evaluated in order, for possible efficiency gain)

  //  evLoader.AddCut("NChan>20");  // example of additional cut
  evLoader.AddCut(testCut);



  //
  // SET NAMES OF VARIABLES TO BE USED FROM ROOT FILES
  //

  TString recoZenRadName = "mZr"; // needed for eProb below
  evLoader.SetName_recoZenith_rad(recoZenRadName);
  // !!!  evLoader.SetName_recoZenith_rad("g32Zr");
  evLoader.SetName_recoAzimuth_rad("mAr");
  //evLoader.SetName_sigmaDeg("mpfSigmaDeg");
  gROOT->ProcessLine(".L mpfSigmaDegRescaled.C+");
  evLoader.SetName_sigmaDeg(" mpfSigmaDegRescaled(mpfSigmaDeg,mmueEn)");
  evLoader.SetName_runID("RunID");
  //evLoader.SetName_runID("Run");
  evLoader.SetName_eventID("EventID");
  // EXAMPLE: How to load different vars for source and bkg, if necessary
  //  dsBkg->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg")
  //  dsSrc->tree->SetAlias("myRecoSigmaDeg","pf32SigmaDeg*1.4");
  //  evLoader.SetName_sigmaDeg("myRecoSigmaDeg");


  // The choice here has consequences for how to fill eProb, below
  if (0) { TString energyVar = "NChan"; }
  if (1) { TString energyVar = "log10(mueEn)"; }

  evLoader.SetName_energyValue(energyVar);



  //
  // LOAD EVENTS, SET UP ENERGY PDFs
  // 

  cout << "\nLoading Background Events -- " << baseTree->GetTitle() << endl;
  cout << "Using cut: " << evLoader.GetCuts()->GetTitle() << endl;
  vector<I3Event> baseEvents;
  evLoader.LoadBkgEvents(baseEvents);

  cout << "\nFilling Energy PDFs:\n";
  // !!!
  if (1) {
    ZenithEnergyProb eProb;
    eProb.SetSourceZenWidthDeg(sourceZenWidthDeg);
    eProb.SetName_recoZenith_rad(recoZenRadName);
    {
      // default: 1 zenith band for whole sky (equiv. to SimpleEnergyProb)
      if (0) { }
      if (0) {
        vector<double> zenMinDegVect;
        zenMinDegVect.push_back(0.);
        zenMinDegVect.push_back(105.);
        zenMinDegVect.push_back(180.);
        eProb.SetZenithBandsDeg(zenMinDegVect);
      }
      if (1) {
        vector<double> zenMinDegVect;

        // Fill vector with bin edges that match "CutDMS" bin edges
        // 0 - 90  zen added
        double tempBot=1;
        for(int i=0; i<21;i++){
          //cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
          zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
          tempBot-=0.05;
        }

        // A reasonable split for the upgoing range
        zenMinDegVect.push_back(105.00); //? Follow Chad's example
        zenMinDegVect.push_back(180.);
        eProb.SetZenithBandsDeg(zenMinDegVect);
      }
    }
  }
  if (0) {
    SimpleEnergyProb eProb;
  }
  // !!!

  eProb.SetLoadModeNew(true); // true is now faster

  if (1) { } // new default is to constrain signal inside histogram
  if (0) { eProb.SetConstrainSignal(false); }  // the old way

  if (energyVar == "NChan") {
    int nBackFill = 20; // don't backfill previous bins
    eProb.SetEnergyGammaRangeAndBackFill(50,0.,1200., 30,1.,4., nBackFill);
  }
  if (energyVar == "log10(mueEn)") {
    int nBackFill = 35; // don't backfill previous bins
    eProb.SetEnergyGammaRangeAndBackFill(40,2.,9., 30,1.,4., nBackFill);
  }

  eProb.SetTableBkg(baseEvents);
  TStopwatch ts;
  eProb.SetTableGamma(srcTree, evLoader.GetCuts(), energyVar);
  ts.Print();

  //  gROOT->ProcessLine(".x example_scripts/monitor_eProb.C");



  // 
  // STORE EVERYTHING IN psData
  //

  I3Analysis psData;  
  psData.SetDecMapParameters(180,1); // 1. deg map bins and smoothing
  // This seems to be okay for ~ 5000 events or more.
  psData.SetBaseEvents(baseEvents);
  psData.SetEnergyProb(eProb);

  // Signal and mySignalPtr can be reloaded and reset anytime
  /*
  double spectralIndex = -2;  // global... can be changed later
  gROOT->ProcessLine(".x macro_loadSignal.C");
  psData.SetSource(*mySignalPtr);
  */

  // for backward compatibility, until other scripts are upgraded:
  //  EventLoader &ns = evLoader;
  //  vector<I3Event> &bkgEvents = baseEvents;

}
