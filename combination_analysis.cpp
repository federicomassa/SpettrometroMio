#include "combinatorial.h"
#include "B_event_approx.h"
#include <TFile.h>
#include <TNtupleD.h>
#include <TH1F.h>
#include <iostream>
 
void combination_analysis() {

  TFile* in = new TFile("../Spettrometro_Files/measures.root");
  TFile* out = new TFile("../Spettrometro_Files/measures_config.root", "recreate");
  TNtupleD* nt = (TNtupleD*) in->Get("measures");
  // TNtupleD* nt2 = new TNtupleD("measures_config", "Best combination", "best_configx:best_configy:chi2_x:chi2_y:a1_x:b1_x:a2_x:b2_x:a1_y:b1_y:a2_y:b2_y");
  TNtupleD* nt2 = new TNtupleD("measures_config", "Best combination", "best_config:K_z:p_plus:p_min:z_reco:p_plus_reco:p_min_reco");
  TH1F* config_hist = new TH1F("config_hist", "Configuraton Histogram", 16, 0, 16);
  Double_t* final_meas = new Double_t[16];
  Double_t* final_pars = new Double_t[2];
  Double_t* args = new Double_t[nt->GetNbranches()];
  Double_t* x1_init = new Double_t[4];
  Double_t* x2_init = new Double_t[4];
  Double_t* y1_init = new Double_t[4];
  Double_t* y2_init = new Double_t[4];
  Double_t* x1_final = new Double_t[4];
  Double_t* x2_final = new Double_t[4];
  Double_t* y1_final = new Double_t[4];
  Double_t* y2_final = new Double_t[4];
  Double_t* z = new Double_t[4];
  Double_t* best_config = new Double_t[7];
  combinatorial* comb;

  z[0] = 50;
  z[1] = 60;
  z[2] = z[1] + B_event_approx::L;
  z[3] = z[2] + B_event_approx::Delta_z;

  for (UInt_t i = 0; i < 5000; i++) {
  // for (UInt_t i = 0; i < 100; i++) {
    if (i % 5000 == 0) std::cout << Double_t(i)/Double_t(nt->GetEntries())*100 << "% completed..." << std::endl;

    comb = new combinatorial;
    
    nt->GetEntry(i);

    args = nt->GetArgs();
    
    for (UInt_t i = 0; i < 4; i++) {
      x1_init[i] = args[20 + 2*i];
      x2_init[i] = args[28 + 2*i];
      y1_init[i] = args[21 + 2*i];
      y2_init[i] = args[29 + 2*i];
    }

    
    
    comb->SetPoints(4, z, x1_init, y1_init, x2_init, y2_init);
    comb->GetCombination(final_meas, final_pars);
    best_config[0] = Double_t(comb->GetBestConfig());
    best_config[1] = args[0];
    best_config[2] = args[2];
    best_config[3] = args[3];
    best_config[4] = final_pars[0];
    best_config[5] = final_pars[1];
    best_config[6] = final_pars[2];
    // best_config[2] = comb->GetBestChi2();
    // comb->GetBest_ab(best_config[4], best_config[5], best_config[6], best_config[7]);
    // best_config[1] = Double_t(comb->GetBestConfig());
    // best_config[3] = comb->GetBestChi2();
    // comb->GetBest_ab(best_config[8], best_config[9], best_config[10], best_config[11]);

    nt2->Fill(best_config);
    delete comb;
    
  }
  nt2->Write();
  in->Close();
  out->Close();

}
  
