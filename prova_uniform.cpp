#include <TRandom3.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>

void prova_uniform() {

  TRandom3 rndgen;
  const int imax = 1E8;
  const double pi = 4*TMath::ATan(1);
  TH1F* prova_hist = new TH1F("P", "Prova", 1000,1,0);
  TFile* out = new TFile("prova_uniform.root","RECREATE");

  for (int i = 0; i < imax; i++) {
    prova_hist->Fill(rndgen.Uniform(2*pi));
  }
  
  prova_hist->Write();
  out->Close();
  
}
