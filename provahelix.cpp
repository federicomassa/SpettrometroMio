#include "helix.h"

#include <TGraph2D.h>
#include <iostream>

void provahelix() {

  Double_t theta = 30*TMath::Pi()/180;
  Double_t phi = 0*TMath::Pi()/180;
  Double_t init_p[3] = {0,0,0};
  Double_t init_sp[3] = {TMath::Sin(theta)*TMath::Cos(phi), TMath::Sin(theta)*TMath::Sin(phi), TMath::Cos(theta)};
  Double_t gamma = 360;
  const int npoints = 10000;
  Double_t r[3];

  TGraph2D* g = new TGraph2D(npoints);
  g->SetTitle("Helix;x;y;z");

  helix h1;
  h1.SetHelix(init_p, init_sp, -1.0, gamma);

  cout << h1.GetRadius() << endl;

  cout << h1.GetGamma() << endl;
  cout << h1.GetCharge() << endl;
  cout << h1.GetOmega() << endl;
  
  for (int i = 0; i < npoints; i++) {
    h1.GetPoint(Double_t(i)*TMath::Power(10,-9),r);
    g->SetPoint(i+1, r[0], r[1], r[2]);
  }

  g->Draw("P");

}
