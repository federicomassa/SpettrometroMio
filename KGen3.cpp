#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>

#include <iostream>

void KGen3() {

  Double_t K_Mass = 0.497614;
  Double_t Pi_ch_Mass = 0.13957018;
  Double_t Pi_neutr_Mass = 0.1349766;
  Double_t K_Energy = 100; //GeV
  Int_t imax = 1E6;
  TLorentzVector beam = {0,0,TMath::Sqrt(K_Energy*K_Energy-K_Mass*K_Mass),K_Energy};

  Double_t mass[3] = {Pi_ch_Mass, Pi_ch_Mass, Pi_neutr_Mass};

  TFile* root_out = new TFile("phacespace_3.root","recreate");
  TH1F* p_hist = new TH1F("p_hist", "Probability of event", 50, 1,0);

  TLorentzVector* pPip = new TLorentzVector;
  TLorentzVector* pPim = new TLorentzVector;
  TLorentzVector* pPi0 = new TLorentzVector;

  TGenPhaseSpace event;
  event.SetDecay(beam,3,mass);

  Double_t weight;
  Double_t weight_max = 0.48; //previosly measured from weight hist
  Double_t prob;
  for (Int_t i = 0; i < imax; i++) {
    
    weight = event.Generate();
    prob = weight/weight_max;
    pPip = event.GetDecay(0);
    pPim = event.GetDecay(1);
    pPi0 = event.GetDecay(2);
  }

  std::cout << weight_sum << std::endl;
  std::cout << event.GetWtMax() << std::endl;
  std::cout << weight_max << std::endl;
  p_hist->DrawCopy();
  root_out->Close();

}
