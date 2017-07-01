
{
  bool optBatch;   // if set to true earlier, will skip drawing plots

  // SET RANGE
  double raLow  = 280;
  double raHigh = 285.;
  double decLow  = -2.;
  double decHigh = 3.;


  // SET COARSE BINNING
  int nBinsCoarseRa = (raHigh-raLow)*2;
  int nBinsCoarseDec = (decHigh-decLow)*2;

  // SET FINE BINNING
  int nBinsFineRa = (raHigh-raLow)*10;
  int nBinsFineDec = (decHigh-decLow)*10;

  // SET -LOG_P VALUE FOR FINE-BINNED FOLLOW-UP
  double resultThreshold = 0.0;



  TH2D *hAllSkyCoarse = new TH2D("hAllSkyCoarse","hAllSkyCoarse",
			   nBinsCoarseRa,raLow,raHigh,
			   nBinsCoarseDec,decLow,decHigh);
  cout << "Coarse Grid: RA Bins: " << nBinsCoarseRa;
  cout << "   Dec Bins: " << nBinsCoarseDec << endl;


  TH2D *hAllSkyFine = new TH2D("hAllSkyFine","hAllSkyFine",
			       nBinsFineRa,raLow,raHigh,
			       nBinsFineDec,decLow,decHigh);
  cout << "Fine Grid:   RA Bins: " << nBinsFineRa;
  cout << "   Dec Bins: " << nBinsFineDec << endl;

  cout << "\nResult Threshold for fine-grid follow-up: " << resultThreshold;
  cout << "\n";


  double decDegMax;
  double raDegMax;
  double resultMax;
  double nsBest;
  double gammaBest;
  double llhBest;

  //
  // COARSE GRID
  //

  TStopwatch ts;
  resultMax = 0.;
  for (int iDec=1; iDec<=nBinsCoarseDec; ++iDec) {
    double decDeg = hAllSkyCoarse->GetYaxis()->GetBinCenter(iDec);
    if (!optBatch) { cout << decDeg << "  " << flush; }

    for (int iRa=1; iRa<=nBinsCoarseRa; ++iRa) {
      double raDeg = hAllSkyCoarse->GetXaxis()->GetBinCenter(iRa);
      
      EquatorialDeg searchLocation(raDeg,decDeg);
      llhFn.SetAnalysis(psData,searchLocation);
      llhFn.MaximizeLlh();

      double result = -log10(llhFn.GetEstProb());
      hAllSkyCoarse->SetBinContent(iRa,iDec, result);

      if (result>resultMax) { 
	resultMax = result;
	decDegMax = decDeg;
	raDegMax = raDeg;
	nsBest = llhFn.GetPar(0);
	gammaBest = llhFn.GetPar(1);
	llhBest = llhFn.Get_logLambdaBest();
      }
    }
  }
  cout << endl;
  ts.Print();

  if (!optBatch) {
    gROOT->SetStyle("Plain");
    CreatePalette(5);
    gStyle->SetOptStat(0);
  
    TCanvas *canAllSky = new TCanvas("canAllSky","canAllSky",800,400);
    hAllSkyCoarse->Draw("colz");
    canAllSky->Update();
  }  

  cout << "Coarse Grid Hottest Spot:\n";
  cout << "   Ra: " << raDegMax << " , Dec: " << decDegMax;
  cout << "   Result Max:  " << resultMax << endl;
  cout << "   llhBest =  " << llhBest << endl;
  cout << "   ns = " << nsBest << "   gamma = " << gammaBest << "\n\n";


  //
  // FINE GRID
  //

  TStopwatch ts;

  // keep same resultMax from above:
  // in case no coarse bin is over threshold, the max results from above
  // still apply

  for (int iDec=1; iDec<=nBinsFineDec; ++iDec) {
    double decDeg = hAllSkyFine->GetYaxis()->GetBinCenter(iDec);
    if (!optBatch) { cout << decDeg << "  " << flush; }

    for (int iRa=1; iRa<=nBinsFineRa; ++iRa) {
      double raDeg = hAllSkyFine->GetXaxis()->GetBinCenter(iRa);

      double result;

      // Check if follow-up required
      double coarseResult = 
	hAllSkyCoarse->GetBinContent( hAllSkyCoarse->FindBin(raDeg,decDeg) );

      if (coarseResult > resultThreshold) {
	EquatorialDeg searchLocation(raDeg,decDeg);
	llhFn.SetAnalysis(psData,searchLocation);
	llhFn.MaximizeLlh();
	result = -log10(llhFn.GetEstProb());

	if (result>resultMax) { 
	  resultMax = result;
	  decDegMax = decDeg;
	  raDegMax = raDeg;
	  nsBest = llhFn.GetPar(0);
	  gammaBest = llhFn.GetPar(1);
	  llhBest = llhFn.Get_logLambdaBest();
	}
      } else {
	result = coarseResult;
      }

      hAllSkyFine->SetBinContent(iRa, iDec, result);

    }
  }

  cout << endl;
  ts.Print();

  if (!optBatch) {
    TCanvas *canAllSkyFine = 
      new TCanvas("canAllSkyFine","canAllSkyFine",800,400);
    hAllSkyFine->Draw("colz");
    canAllSkyFine->Update();
  }

  cout << "Fine Grid Hottest Spot:\n";
  cout << "Ra: " << raDegMax << " , Dec: " << decDegMax;
  cout << "   Result Max:  " << resultMax << endl;
  cout << "   llhBest =  " << llhBest << endl;
  cout << "   ns = " << nsBest << "   gamma = " << gammaBest << endl;

}
