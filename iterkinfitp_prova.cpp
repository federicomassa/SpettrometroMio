#include "IterKinFitP.h"
#include "constr_prova.h"

#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>

void iterkinfitp_prova() {

  TRandom3 rndgen;
  IterKinFitP iter;
  TFile* out = new TFile("iterkinfitp_prova.root", "recreate");
  Double_t theta = 0.005;
  Double_t phi = 45*TMath::Pi()/180;
  Double_t z_true = 15;
  Double_t* err = new Double_t[4];

  err[0] = err[1] = err[2] = err[3] = 0.001;

  Double_t* init_meas_t = new Double_t[4];
  Double_t* init_pars_t = new Double_t;

  init_meas_t[0] = (z1 - z_true)*TMath::Tan(theta)*TMath::Cos(phi);
  init_meas_t[1] = (z1 - z_true)*TMath::Tan(theta)*TMath::Sin(phi);
  init_meas_t[2] = (z2 - z_true)*TMath::Tan(theta)*TMath::Cos(phi);
  init_meas_t[3] = (z2 - z_true)*TMath::Tan(theta)*TMath::Sin(phi);

  *init_pars_t = z_true;


  Double_t* init_meas = new Double_t[4];
  Double_t* init_pars = new Double_t;

  init_meas[0] = rndgen.Gaus(init_meas_t[0], err[0]);
  init_meas[1] = rndgen.Gaus(init_meas_t[1], err[1]);
  init_meas[2] = rndgen.Gaus(init_meas_t[2], err[2]);
  init_meas[3] = rndgen.Gaus(init_meas_t[3], err[3]);

  *init_pars = (init_meas[2]*z1 - init_meas[0]*z2)/(init_meas[2]-init_meas[0]);

  std::cout << "TRUE MEAS: " << std::endl;

  for (UInt_t i = 0; i < 4; i++) 
    std::cout << init_meas_t[i] << std::endl;

  std::cout << *init_pars_t << std::endl;
  std::cout << std::endl;

  iter.Initialize(4,2,1, init_meas, init_pars, Constr, Der_Matrix, PDer_Matrix, err);

  TGraph* graph = new TGraph[5];

  Double_t* final_meas = new Double_t[4];
  Double_t* final_pars = new Double_t;

  iter.Minimize(final_meas, final_pars, graph);

  iter.PrintResult();

  for (int i = 0; i < 5; i++) {
    graph[i].Write();}
  out->Close();
  
}
