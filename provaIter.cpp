#include "IterKinFit.h"
#include <TMatrixD.h>
#include <TRandom3.h>
#include <iostream>
#include <iomanip>
#include <TFile.h>

TMatrixD Constraints(Double_t* var) {
  
  TMatrixD Phi(1,1);
  
  Phi(0,0) = var[0] + var[1] + var[2] - 180;
  return Phi;

}

TMatrixD Derivatives(Double_t* var) {

  TMatrixD Der(3,1);
  
  Der(0,0) = 1;
  Der(1,0) = 1;
  Der(2,0) = 1;

  return Der;
}


void provaIter() {
  
  TFile* out = new TFile("provaIter.root","recreate");
  TRandom3 rndgen;
  Double_t err[3] = {0.001, 0.001, 0.001};
  IterKinFit iter;

  Double_t init[3] = {rndgen.Gaus(70,0.001), rndgen.Gaus(80, 0.001), rndgen.Gaus(30,0.001)};
  Double_t final[3];

  iter.Initialize(3,1, init, Constraints, Derivatives, err);
  iter.SetStepParameter(1.0);
  iter.SetThreshold(1E-6);
  TGraph gr[3];
  iter.Minimize(final, gr);

  std::cout << std::fixed << std::setprecision(10) << final[0] << '\n' << final[1] << '\n' << final[2] << std::endl;

  std::cout << "iter_number: " << iter.GetIterationNumber() << std::endl;

  gr[0].Write();
  gr[1].Write();
  gr[2].Write();

  out->Close();
  

}
