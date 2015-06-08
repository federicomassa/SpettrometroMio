#include <TH1F.h>
#include <TFile.h>

void provaroot() {

  TFile* out = new TFile("provaroot.root", "recreate");
  TH1F* h1 = new TH1F("h1", "HIST", 1000,-10,10);
  TH1F* h2 = new TH1F("h2", "HIST", 1000,-10,10);

  h1->FillRandom("gaus",10000);
  h2->FillRandom("gaus",10000);
  h1->Write();
  h2->Write();
  out->Close();
  

}
