#include "IterKinFitP.h"
#include "const_double2_1det.h"
#include "InitializeZPP_1det.h"
#include "MergeNtuples.h"
#include <TNtupleD.h>
#include <TFile.h>
#include <iostream>
#include <TMath.h>

void p_reco() {


  const char* nt_iter_list = "chi2:z_r:p_plus_r:p_min_r:x1_plus_r:y1_plus_r:x2_plus_r:y2_plus_r:x3_plus_r:y3_plus_r:x1_min_r:y1_min_r:x2_min_r:y2_min_r:x3_min_r:y3_min_r:W_r";

  TFile* root_in = new TFile("../Spettrometro_Files/measures.root");
  TFile* root_out = new TFile("../Spettrometro_Files/measures_final_double2_1det.root", "recreate");
  TNtupleD* nt_gen = (TNtupleD*) root_in->Get("measures");
  TNtupleD* nt_iter = new TNtupleD("nt_iter", "Iterative KinFit", nt_iter_list);
  TNtupleD* nt_final = new TNtupleD;
  IterKinFitP* iter;
  Double_t* args_gen = new Double_t[36];
  Double_t* init_meas = new Double_t[12];
  Double_t* init_pars = new Double_t[3];
  Double_t* final_meas = new Double_t[12];
  Double_t* final_pars = new Double_t[3];
  Double_t* err = new Double_t[12];
  Double_t chi2;
  Double_t* args_iter = new Double_t[17];
  for (UInt_t i = 0; i < 12; i++) {
    err[i] = 0.001;
  }
  

  for (UInt_t i = 0; i < nt_gen->GetEntries(); i++) {
  // for (UInt_t i = 0; i < 2000; i++) {
    if (i % TMath::Nint(double(nt_gen->GetEntries())/100*1) == 0) std::cout << TMath::Nint(double(i)/double(nt_gen->GetEntries())*100) << "% completed..." << std::endl;

    nt_gen->GetEntry(i);
    args_gen = nt_gen->GetArgs();

    for (UInt_t j = 0; j < 6; j++) {
      init_meas[j] = args_gen[20 + j];
    }

    for (UInt_t j = 0; j < 6; j++) {
      init_meas[j+6] = args_gen[28 + j];
    }

    InitializeZPP_1det(init_meas, init_pars);

    // std::cout << init_pars[0] << '\t' << init_pars[1] << '\t' << init_pars[2] << std::endl;

    // std::cout << "////////////////////////" << std::endl;
    // std::cout << i << std::endl;

    iter = new IterKinFitP;
    iter->Initialize(12,11,3,init_meas,init_pars,Constr, Der, PDer, err);
    // iter->SetStepParameter(0.5);
    // iter->SetThreshold(1E-5);
    // iter->SetMaxIterationNumber(2000);
    iter->Minimize(final_meas, final_pars);
    chi2 = iter->GetChi2();

    // std::cout << "CHI2: " << chi2 << std::endl;

    if (iter->IsOverMax()) std::cout << "WARNING: Max Iteration Number exceeded" << std::endl;

    args_iter[0] = chi2;
    args_iter[1] = final_pars[0];
    args_iter[2] = 1.0/final_pars[1];
    args_iter[3] = -1.0/final_pars[2]; // negative charge

    for (UInt_t j = 4; j < 16; j++) {
      args_iter[j] = final_meas[j-4];
    }

    args_iter[16] = Invariant_Mass(final_meas, final_pars);

    nt_iter->Fill(args_iter);

  }


  nt_final =  MergeNtuples(nt_gen, nt_iter, "measures_final_double2_1det", "Complete Ntuple");

  nt_final->Write();
  root_in->Close();
  root_out->Close();

}
