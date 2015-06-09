#include "IterKinFit.h"
#include "momentum_constraint.h"

#include <TMath.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>

void provamomentum() {

  IterKinFit iter;
  TRandom3 rnd;
  TFile* out = new TFile("provamomentum.root", "recreate");
  Double_t err[6] = {0.001, 0.001, 1, 0.001, 0.001, 1};
  Double_t true_meas[6] = {-0.0309593, 0.148942 , 23.1974, 0.0309593, -0.148942, 81.7973};



  Double_t init[6] = {rnd.Gaus(true_meas[0],err[0]), rnd.Gaus(true_meas[1],err[1]), rnd.Gaus(true_meas[2],err[2]), rnd.Gaus(true_meas[3],err[3]), rnd.Gaus(true_meas[4],err[4]), rnd.Gaus(true_meas[5],err[5])};

  std::cout << "Initial error in W: " << (Momentum_Constr_Vector(init))(2,0) << std::endl;


  iter.Initialize(6,3,init, Momentum_Constr_Vector, Momentum_Constr_Derivative, err);

  iter.SetThreshold(1E-8);
  
  Double_t final[6];
  TGraph* graph = new TGraph[6];
  iter.Minimize(final, graph);
  iter.PrintResult();
  
  for (UInt_t i = 0; i < 6; i++) {
    graph[i].Write();
  }

  std::cout << "Final error in W: " << (Momentum_Constr_Vector(final))(2,0) << std::endl;

  out->Close();

}
