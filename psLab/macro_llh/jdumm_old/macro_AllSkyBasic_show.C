
{
  gROOT->SetStyle("Plain");
  CreatePalette(5);

  // Large Size
  can = new TCanvas("AllSky","AllSky",20,20,1500,600);

  // Medium Size
  //    can = new TCanvas("AllSky","AllSky",20,20,1000,400);

  can->cd();
    
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.075); 
    
    TH2D *h = hAllSkyFine;
    h->SetTitle("");  
    TH2D* hAllSkyFine_surface = dynamic_cast<TH2D*> h->DrawCopy("surf2z");
    //  TH2D* hAllSkyFine_surface = dynamic_cast<TH2D*> h->DrawCopy("colz");
    gPad->SetTheta(270.);
    gPad->SetPhi(179.99);
    
    hAllSkyFine_surface->SetStats(kFALSE);
    TAxis *axis;
    axis = hAllSkyFine_surface->GetXaxis();
    axis->SetNdivisions(-20306);
    axis->SetLabelSize(.05);
    axis->SetLabelOffset(0.02);
    axis->SetTitle("right ascension");
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(1.2);
    axis->CenterTitle();
    
    axis = hAllSkyFine_surface->GetYaxis();
    axis->SetLabelSize(.05);
    axis->SetTickLength(0.01);
    axis->SetTitle("declination");
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(0.6);
    axis->CenterTitle();
    can->Update();

    TPaletteAxis *palette = dynamic_cast<TPaletteAxis*> 
      (hAllSkyFine_surface->GetListOfFunctions()->FindObject("palette"));

    palette->GetAxis()->SetTitle("-log_{10} p");
    palette->GetAxis()->CenterTitle();
    palette->GetAxis()->SetTitleOffset(0.4);
    palette->GetAxis()->SetTitleSize(0.05);
    palette->Draw();
    can->Update();


  bool opt_Events = 1;
  if (opt_Events) {
    TGraph *eventGraph = Graph_I3Events_Equatorial(psData->GetEventList());
    if (gPad->GetView()) { // 3D
      TPolyMarker3D *eventPoly = GraphToPolyMarker3D(eventGraph,6);
      eventPoly->Draw();
    } else {
      eventGraph->SetMarkerStyle(6);
      eventGraph->Draw("P");
    }
  }

  bool opt_Sources = 1;
  if (opt_Sources) {
    if (srcGraph) {
      if (gPad->GetView()) { // 3D
        TPolyMarker3D *srcPoly = GraphToPolyMarker3D(srcGraph,2);
        srcPoly->SetMarkerColor(2);
        srcPoly->Draw();
      } else {
        srcGraph->SetMarkerStyle(2);
        srcGraph->SetMarkerColor(2);
        srcGraph->Draw("P");
      }
    } else {
      cout << "No srcGraph\n";
    }
  }

  bool opt_GalPlane = 1;
  if (opt_GalPlane) {
    TGraph *galPlane = new TGraph("graphics/GalPlane3600.coords");
    if (gPad->GetView()) { // 3D
      TPolyMarker3D *galPoly = GraphToPolyMarker3D(galPlane,1);
      galPoly->Draw();
    } else {
      galPlane->SetMarkerStyle(1);
      galPlane->Draw("P");
    }
  }

}
