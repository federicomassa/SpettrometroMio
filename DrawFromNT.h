//especially written for ntuple with 36 branches containing all the event information

#ifndef DRAWFROMNT_H
#define DRAWFROMNT_H

#ifndef ZARRAY
#define ZARRAY
#include "B_event_approx.h"
Double_t z1 = 50;
Double_t z2 = 60;
Double_t z3 = z2 + B_event_approx::L;
Double_t z4 = z3 + B_event_approx::Delta_z;
#endif

#include <TNtupleD.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <iostream>
#include <TCanvas.h>
#include "InitializeZPP_1det.h"
#include "const_double2_1det.h"
#include "IterKinFitP.h"


void DrawFromNt(TNtupleD* nt, Int_t event_number, Int_t detector_number) {

  nt->GetEntry(event_number);
  Double_t* args = new Double_t[36];
  args = nt->GetArgs();
  
  Double_t* x_plus = new Double_t[detector_number];
  Double_t* y_plus = new Double_t[detector_number];
  Double_t* x_min = new Double_t[detector_number];
  Double_t* y_min = new Double_t[detector_number];
  Double_t* z = new Double_t[detector_number];
 

  TGraph* track1XZ = new TGraph;
  TGraph* track2XZ = new TGraph;
  TGraph* track1YZ = new TGraph;
  TGraph* track2YZ = new TGraph;

  TMultiGraph* mgXZ = new TMultiGraph();
  TMultiGraph* mgYZ = new TMultiGraph();
  
  if (detector_number == 3) {
    z[0] = z1;
    z[1] = z2;
    z[2] = z3;
  }
  else if (detector_number == 4) {
    z[0] = z1;
    z[1] = z2;
    z[2] = z3;
    z[3] = z4;
  }
  else {std::cout << "ERROR in DrawFromNt, how many detectors?" << std::endl;}
  
  for (Int_t i = 0; i < detector_number; i++) {
    x_plus[i] = args[20 + 2*i];
    y_plus[i] = args[21 + 2*i];   
    x_min[i] = args[28 + 2*i];
    y_min[i] = args[29 + 2*i];   
    track1XZ->SetPoint(i, x_plus[i], z[i]);
    track1YZ->SetPoint(i,y_plus[i], z[i]);
    track2XZ->SetPoint(i, x_min[i], z[i]);
    track2YZ->SetPoint(i,y_min[i], z[i]);
  }

  track1XZ->SetMarkerStyle(kFullCircle);
  track1XZ->SetMarkerColor(kRed);
  
  track2XZ->SetMarkerStyle(kFullCircle);
  track2XZ->SetMarkerColor(kBlue);

  track1YZ->SetMarkerStyle(kFullCircle);
  track1YZ->SetMarkerColor(kRed);
  
  track2YZ->SetMarkerStyle(kFullCircle);
  track2YZ->SetMarkerColor(kBlue);
  
  mgXZ->Add(track1XZ);
  mgXZ->Add(track2XZ);
  mgYZ->Add(track1YZ);
  mgYZ->Add(track2YZ);
  
  TCanvas* c1 = new TCanvas("XZ", "XZ");
  mgXZ->Draw("AP");
  TCanvas* c2 = new TCanvas("YZ", "YZ");
  mgYZ->Draw("AP");

  IterKinFitP* iter = new IterKinFitP;
  Double_t* init_meas = new Double_t[4*detector_number];
  Double_t* init_pars = new Double_t[3];
  Double_t* err = new Double_t[4*detector_number];
  Double_t* final_meas = new Double_t[4*detector_number];
  Double_t* final_pars = new Double_t[3];

  for (Int_t i = 0; i < 2*detector_number; i++) {
    init_meas[i] = args[20+i];
    err[i] = 0.001;
  }
  for (Int_t i = 0; i < 2*detector_number; i++) {
    init_meas[i+2*detector_number] = args[28 +i];
    err[i+2*detector_number] = 0.001;
  }
 

  InitializeZPP_1det(init_meas, init_pars);
  iter->Initialize(4*detector_number, 15 - 4*(4-detector_number), 3, init_meas, init_pars, Constr, Der, PDer, err);
  iter->Minimize(final_meas, final_pars);
  iter->PrintResult();
}

#endif
