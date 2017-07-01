void macro_AllSkyNew_save(char *filename) {
  TFile *f = new TFile(filename,"new");

  hAllSkyCoarse->Write();
  hAllSkyFine->Write();

  f->Close();
}

