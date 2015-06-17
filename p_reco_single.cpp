#include "IterKinFitP.h"
#include "const_single.h"
#include "InitializeZP.h"
#include "MergeNtuples.h"
#include <TNtupleD.h>
#include <TFile.h>
#include <iostream>
#include <TMath.h>

void p_reco() {


  const char* nt_iter_plus_list = "chi2_plus:z_plus_r:p_plus_r:x1_plus_r:y1_plus_r:x2_plus_r:y2_plus_r:x3_plus_r:y3_plus_r:x4_plus_r:y4_plus_r";
  const char* nt_iter_min_list = "chi2_min:z_min_r:p_min_r:x1_min_r:y1_min_r:x2_min_r:y2_min_r:x3_min_r:y3_min_r:x4_min_r:y4_min_r";

  TFile* root_in = new TFile("../Spettrometro_Files/measures.root");
  TFile* root_out = new TFile("../Spettrometro_Files/measures_final_single.root", "recreate");
  TNtupleD* nt_gen = (TNtupleD*) root_in->Get("measures");
  TNtupleD* nt_iter_plus = new TNtupleD("nt_iter_plus", "Iterative KinFit", nt_iter_plus_list);
  TNtupleD* nt_iter_min = new TNtupleD("nt_iter_min", "Iterative KinFit", nt_iter_min_list);
  TNtupleD* nt_iter = new TNtupleD;;
  TNtupleD* nt_final = new TNtupleD;
  IterKinFitP* iter;
  Double_t* args_gen = new Double_t[36];
  Double_t* init_meas = new Double_t[8];
  Double_t* init_pars = new Double_t[2];
  Double_t* final_meas = new Double_t[8];
  Double_t* final_pars = new Double_t[2];
  Double_t* err = new Double_t[8];
  Double_t chi2;
  Double_t* args_plus = new Double_t[11];
  Double_t* args_min = new Double_t[11];
  Double_t* args_final = new Double_t[22];
  for (UInt_t i = 0; i < 8; i++) {
    err[i] = 0.001;
  }
  

  for (UInt_t i = 0; i < nt_gen->GetEntries(); i++) {

    if (i % TMath::Nint(double(nt_gen->GetEntries())/100*5) == 0) std::cout << TMath::Nint(double(i)/double(nt_gen->GetEntries())*100) << "% completed..." << std::endl;

    nt_gen->GetEntry(i);
    args_gen = nt_gen->GetArgs();

    for (UInt_t j = 0; j < 8; j++) {
      init_meas[j] = args_gen[20 + j];
    }

    InitializeZP(init_meas, init_pars);

    iter = new IterKinFitP;
    iter->Initialize(8,6,2,init_meas,init_pars,Constr, Der, PDer, err);
    iter->Minimize(final_meas, final_pars);
    chi2 = iter->GetChi2();

    args_plus[0] = chi2;
    args_plus[1] = final_pars[0];
    args_plus[2] = 1.0/final_pars[1];

    for (UInt_t j = 3; j < 11; j++) {
      args_plus[j] = final_meas[j-3];
    }

    delete iter;

    for (UInt_t j = 0; j < 8; j++) {
      init_meas[j] = args_gen[28 + j];
    }
    
    InitializeZP(init_meas, init_pars);

    iter = new IterKinFitP;
    iter->Initialize(8,6,2,init_meas,init_pars,Constr, Der, PDer, err);
    iter->Minimize(final_meas, final_pars);
    chi2 = iter->GetChi2();

    args_min[0] = chi2;
    args_min[1] = final_pars[0];
    args_min[2] = 1.0/final_pars[1];

    for (UInt_t j = 3; j < 11; j++) {
      args_min[j] = final_meas[j-3];
    }

    nt_iter_plus->Fill(args_plus);
    nt_iter_min->Fill(args_min);

  }

  nt_iter = MergeNtuples(nt_iter_plus, nt_iter_min, "measures_iter", "Iter Ntuple");
  nt_final =  MergeNtuples(nt_gen, nt_iter, "measures_final", "Complete Ntuple");

  nt_final->Write();
  root_in->Close();
  root_out->Close();

}
