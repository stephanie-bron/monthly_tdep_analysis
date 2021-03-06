
{
  // FOR FELDMAN COUSINS
  //  gROOT->ProcessLine(".L FluxTools.C");

  gROOT->ProcessLine(".x macro_loadPointSources.C");

  
  bool batchPointSourceBasic;
  bool OPT_FCLIMITS;
  bool OPT_TEX;
  bool OPT_UPDATE;

  if (batchPointSourceBasic) {
    OPT_FCLIMITS = false;
    OPT_TEX = false;
    OPT_UPDATE = false;
  }
  // if batch not set to true somewhere else, then you can these manually here
  else {
    bool OPT_FCLIMITS = false;
    bool OPT_TEX = false;
    bool OPT_UPDATE = true;
  }


  /*  FC STUFF
  char* fcDirectory = "/net/user/cfinley/npx/sigUpgrade_Em2_Nch/results/";
  double sigUnc = 0.15;  // systematic unc. on signal
  double spectralIndex = -2;  // necessary if you are converting fluxes
  double old_EnUnits = 1e9;
  double new_EnUnits = 1e12;
  double old_EnReference = 1e9;
  double new_EnReference = 1e12;
  */

  
  TString texOutputString;
  TString plainOutputString;

  double bestPMin = 1.;
  int bestSource = 0;
      
  for (int j=0; j<nPsList; ++j) {
    EquatorialDeg searchLocation = psEqCoordArray[j];
    llhFn.SetAnalysis(psData,searchLocation);
    llhFn.MaximizeLlh();
    
    double pvalue = llhFn.GetEstProb();
    if (pvalue < bestPMin) { 
      bestPMin = pvalue; 
      bestSource = j;
    }
    double distance = 2.; // degrees
    double bkgCount = psData.BkgNumberDensity(searchLocation)
      * (distance)*(distance)*TMath::Pi();
    double nSrc = llhFn.GetPar(0);
    
    double fluxUpLimit = 0.;

    if (OPT_FCLIMITS) {
      cout << "FC LIMITS NOT CURRENTLY IMPLEMENTED.\n";
      /*
      SignedSqrtFC signedfc;
      signedfc.SetParams(800,-10.,30.,410, sigUnc);
      SetDecGeneralFC(signedfc,fcDirectory, 50, 30., decDeg);
      double confLev = 0.90;
      signedfc.BuildConfidenceBand(confLev);
      
      double maxLlh = llhFn.Get_logLambdaBest();
      double nSrcUpLimit = signedfc.GetUpperLimit(sqrt(2*fabs(maxLlh)));
      
      ns.SetSourceEqDeg( searchLocation );
      ns.LoadSourceEvents( sourceEvents );
      NugenSource mySource;
      mySource.SetSourceParams(flux, spectralIndex, livetime);
      mySource.StoreCandidateEvents(searchLocation, sourceEvents);
      
      fluxUpLimit = mySource.GetFluxForNev(nSrcUpLimit) *
	ScaleFlux(old_EnUnits, new_EnUnits, old_EnReference,
		  new_EnReference, spectralIndex) *
	1e12; // report things in units of TeV
      
      // int evCount = EventsNearSource(*(psData.GetEvents()),
      //	       		     searchLocation, distance).size();
      
      // srcListDec.push_back(decDeg);
      // srcListFluxUpLim.push_back(mySource.GetFluxForNev(nSrcUpLimit));
      */
    }

    TString sName = psNameArray[j];
    sName.ReplaceAll("_"," ");
      
    char buffer[1000];

    if (OPT_TEX) {
      char pvalueChar[20];
      if (pvalue < 0.495) {
	sprintf(pvalueChar,"%.2lg",pvalue);
      } else {
	sprintf(pvalueChar,"--");
      }
      sprintf(buffer,
	      "%20s & %6.2f & %5.2f & %5.2f & %s & %3.1f & %3.1f\\\\\n",
	      dynamic_cast<char*>(sName),
	      searchLocation.GetRa(), searchLocation.GetDec(),
	      fluxUpLimit,
	      pvalueChar,
	      nSrc,
	      bkgCount);
      texOutputString += buffer;
    }
    
    if (OPT_UPDATE) {
      char pvalueChar[20];
      if (pvalue < 0.495) {
	sprintf(pvalueChar,"%.5lf",pvalue);
      } else {
	sprintf(pvalueChar,"  ---  ");
      }
      sprintf(buffer,
	      "%20s  (%7.3f , %6.3f)  :  %5.2f  %s  %3.1f   %3.1f\n",
	      dynamic_cast<char*>psNameArray[j],
	      searchLocation.GetRa(), searchLocation.GetDec(),
	      fluxUpLimit,
	      pvalueChar,
	      nSrc,
	      bkgCount);
      plainOutputString += buffer;
      printf("%s",buffer);
    }
  }

  if (OPT_TEX) {
    printf("\n\nTEX OUTPUT:\n");
    cout << (char*)texOutputString;
  }
  cout << "Smallest p-value = " <<bestPMin << endl;
  cout << "for Source " << bestSource << " (" << dynamic_cast<char*>(sName) << ")\n";


}
