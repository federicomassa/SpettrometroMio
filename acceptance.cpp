#include <TFile.h>
#include <TGraphErrors.h>
#include <TNtupleD.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <string>
#include <sstream>
#include <iostream>

Double_t R_min = 0.1;
Double_t R_max = 0.2;
Int_t nR = 21;
Double_t DeltaR = (R_max-R_min)/(Double_t(nR)-1);

void acceptance() {
  Int_t accepted;
  Double_t acceptance;
  Int_t generated = 1E5;
  stringstream ss;
  string R_str;
  TFile in("../Spettrometro_Files/measures.root");
  TFile out("../Spettrometro_Files/acceptance.root", "recreate");
  TNtupleD* nt = (TNtupleD*) in.Get("measures");

  TGraphErrors acc;

  for (Int_t i = 0; i < nR; i++) {
    ss << R_min + Double_t(i)*DeltaR;
    ss >> R_str;
    ss.clear();

    accepted = nt->GetEntries(("TMath::Sqrt(x1_plus_t**2+y1_plus_t**2) > " + R_str + "&& TMath::Sqrt(x1_min_t**2+y1_min_t**2) > " + R_str).c_str());
    acceptance = Double_t(accepted)/Double_t(generated);

    acc.SetPoint(i, R_min+Double_t(i)*DeltaR, acceptance);
    acc.SetPointError(i, 0.0, TMath::Sqrt(acceptance*(1-acceptance)/generated));

  }

  acc.Write();
  in.Close();
  out.Close();

}
