#include "B_event.h"
#include "GetThetaPhi.h"
#include "ReadImpr.h"
#include "IterKinFit.h"
#include "momentum_constraint.h"
#include <fstream>
#include <string>
#include <sstream>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TGraph.h>

using namespace TMath;

void p_reco(const char* impr_meas_name = "impr_measures") {

  string path = "../Spettrometro_Files/";
  string dat_ext = ".dat";
  string impr_str;
  string complete_impr_name;
  string str;

  TRandom3 rndgen;

  Double_t P1_plus[3], P1_min[3], P2_plus[3], P2_min[3], P3_plus[3], P3_min[3], P4_plus[3], P4_min[3];
  Double_t P1_plus_t[3], P1_min_t[3], P2_plus_t[3], P2_min_t[3], P3_plus_t[3], P3_min_t[3], P4_plus_t[3], P4_min_t[3];
  Double_t det1 = 50;
  Double_t det2 = 60;
  Int_t ev_no;
  Double_t K_z, K_p, pi_plus_modp, pi_min_modp;
  
  TGraph* graph = new TGraph[6];

  Double_t pi_plus_px, pi_min_px;
  Double_t pi_plus_py, pi_min_py;
  Double_t pi_plus_pz, pi_min_pz;

  Double_t p_reco_init[6], p_reco_final[6], p_reco_err[6], p_true[6];
  p_reco_err[0] = p_reco_err[3] = 0.00365;
  p_reco_err[1] = p_reco_err[4] = 0.00345;
  p_reco_err[2] = p_reco_err[5] = 0.47;

  Double_t gamma_plus, gamma_min;
  Double_t theta_plus_in_t, theta_plus_out_t;
  Double_t phi_plus_in_t;
  Double_t theta_min_in_t, theta_min_out_t;
  Double_t phi_min_in_t;

  Double_t theta_plus_in, theta_plus_out;
  Double_t phi_plus_in, phi_plus_out;
  Double_t theta_min_in, theta_min_out;
  Double_t phi_min_in, phi_min_out;

  Double_t pi_plus_p_reco, pi_min_p_reco;
  Double_t pi_plus_px_reco, pi_min_px_reco;
  Double_t pi_plus_py_reco, pi_min_py_reco;
  Double_t pi_plus_pz_reco, pi_min_pz_reco;

  Double_t pi_plus_p_final, pi_min_p_final;

  Double_t pi_mass = 0.13957018; // [GeV/c^2]
  B_event ev_plus;
  B_event ev_min;
  IterKinFit* iter;

  int count = 0;
  //HISTOGRAMS
  TH1F* p_plus_reco_hist = new TH1F("p_plus_reco", "Pi+ p Reco; True-Meas", 1000, 1, 0);
  TH1F* p_min_reco_hist = new TH1F("p_min_reco", "Pi- p Reco; True-Meas", 1000, 1, 0);
  TH1F* p_plus_final_hist = new TH1F("p_plus_final", "Pi+ p Final; True-Meas", 1000, 1, 0);
  TH1F* p_min_final_hist = new TH1F("p_min_final", "Pi- p Final; True-Meas", 1000, 1, 0);

  TH1F* px_plus_reco_hist = new TH1F("px_plus_reco", "Pi+ px Reco; True-Meas", 1000, 1, 0);
  TH1F* py_plus_reco_hist = new TH1F("py_plus_reco", "Pi+ py Reco; True-Meas", 1000, 1, 0);
  TH1F* pz_plus_reco_hist = new TH1F("pz_plus_reco", "Pi+ pz Reco; True-Meas", 1000, 1, 0);

  TH1F* px_min_reco_hist = new TH1F("px_min_reco", "Pi- px Reco; True-Meas", 1000, 1, 0);
  TH1F* py_min_reco_hist = new TH1F("py_min_reco", "Pi- py Reco; True-Meas", 1000, 1, 0);
  TH1F* pz_min_reco_hist = new TH1F("pz_min_reco", "Pi- pz Reco; True-Meas", 1000, 1, 0);

  TH1F* px_plus_final_hist = new TH1F("px_plus_final", "Pi+ px Final; True-Meas", 1000, 1, 0);
  TH1F* py_plus_final_hist = new TH1F("py_plus_final", "Pi+ py Final; True-Meas", 1000, 1, 0);
  TH1F* pz_plus_final_hist = new TH1F("pz_plus_final", "Pi+ pz Final; True-Meas", 1000, 1, 0);

  TH1F* px_min_final_hist = new TH1F("px_min_final", "Pi- px Final; True-Meas", 1000, 1, 0);
  TH1F* py_min_final_hist = new TH1F("py_min_final", "Pi- py Final; True-Meas", 1000, 1, 0);
  TH1F* pz_min_final_hist = new TH1F("pz_min_final", "Pi- pz Final; True-Meas", 1000, 1, 0);
  TH1F* K_reco_hist = new TH1F("K_reco_hist", "K_true - K_reco; True-Meas", 1000, 1, 0);
  TH1F* K_final_hist = new TH1F("K_final_hist", "K_true - K_reco; True-Meas", 1000, 1, 0);
  TH1F* K_final2_hist = new TH1F("K_final2_hist", "K_true - K_reco; True-Meas", 1000, 1, 0);
  TH1F* K_final3_hist = new TH1F("K_final3_hist", "K_true - K_reco; True-Meas", 1000, 1, 0);


  ////////////

  // Detector
  Double_t sigma = 0.001;


  P1_plus[2] = P1_min[2] = det1;
  P2_plus[2] = P2_min[2] = det2;
  
  stringstream ss;
  ss << impr_meas_name;
  ss >> impr_str;
  ss.clear();
  
  complete_impr_name = path + impr_str + dat_ext;

  ifstream impr_file(complete_impr_name.c_str());
  TFile* root_out = new TFile("../Spettrometro_Files/p_reco.root","recreate");
  
  getline(impr_file,str);
  getline(impr_file,str);
  getline(impr_file,str);
  getline(impr_file,str);
  getline(impr_file,str);
    
  ReadImpr(impr_file, ev_no, K_z, K_p, pi_plus_modp, theta_plus_in_t, phi_plus_in_t, pi_min_modp, theta_min_in_t, phi_min_in_t, P1_plus, P1_min, P2_plus, P2_min);

  do {

    count++;
    if (count % 20000 == 0) cout << count << " events..." << endl;

    P1_plus_t[0] = (det1-K_z)*Tan(theta_plus_in_t)*Cos(phi_plus_in_t);
    P1_plus_t[1] = (det1-K_z)*Tan(theta_plus_in_t)*Sin(phi_plus_in_t);
    P1_plus_t[2] = det1;

    P1_min_t[0] = (det1-K_z)*Tan(theta_min_in_t)*Cos(phi_min_in_t);
    P1_min_t[1] = (det1-K_z)*Tan(theta_min_in_t)*Sin(phi_min_in_t);
    P1_min_t[2] = det1;

    P2_plus_t[0] = (det2-K_z)*Tan(theta_plus_in_t)*Cos(phi_plus_in_t);
    P2_plus_t[1] = (det2-K_z)*Tan(theta_plus_in_t)*Sin(phi_plus_in_t);
    P2_plus_t[2] = det2;

    P2_min_t[0] = (det2-K_z)*Tan(theta_min_in_t)*Cos(phi_min_in_t);
    P2_min_t[1] = (det2-K_z)*Tan(theta_min_in_t)*Sin(phi_min_in_t);
    P2_min_t[2] = det2;

    //True momentum vector
    pi_plus_px = pi_plus_modp*Sin(theta_plus_in_t)*Cos(phi_plus_in_t);
    pi_plus_py = pi_plus_modp*Sin(theta_plus_in_t)*Sin(phi_plus_in_t);
    pi_plus_pz = pi_plus_modp*Cos(theta_plus_in_t);

    pi_min_px = pi_min_modp*Sin(theta_min_in_t)*Cos(phi_min_in_t);
    pi_min_py = pi_min_modp*Sin(theta_min_in_t)*Sin(phi_min_in_t);
    pi_min_pz = pi_min_modp*Cos(theta_min_in_t);

    p_true[0] = pi_plus_px;
    p_true[1] = pi_plus_py;
    p_true[2] = pi_plus_pz;
    p_true[3] = pi_min_px;
    p_true[4] = pi_min_py;
    p_true[5] = pi_min_pz;


    ///////////

    gamma_plus = Sqrt(Power(pi_plus_modp/pi_mass,2)+1);
    ev_plus.SetB_event(P1_plus_t, P2_plus_t, 1.0, gamma_plus);

    gamma_min = Sqrt(Power(pi_min_modp/pi_mass,2)+1);
    ev_min.SetB_event(P1_min_t, P2_min_t, -1.0, gamma_min);


    ev_plus.GetP3(P3_plus_t);
    ev_min.GetP3(P3_min_t);

    ev_plus.GetP4(P4_plus_t);
    ev_min.GetP4(P4_min_t);

    P3_plus[0] = rndgen.Gaus(P3_plus_t[0],sigma);
    P3_plus[1] = rndgen.Gaus(P3_plus_t[1],sigma);
    P3_plus[2] = det2 + B_event::L;

    P3_min[0] = rndgen.Gaus(P3_min_t[0],sigma);
    P3_min[1] = rndgen.Gaus(P3_min_t[1],sigma);
    P3_min[2] = det2 + B_event::L;

    P4_plus[0] = rndgen.Gaus(P4_plus_t[0],sigma);
    P4_plus[1] = rndgen.Gaus(P4_plus_t[1],sigma);
    P4_plus[2] = det2 + B_event::L + B_event::Delta_z;

    P4_min[0] = rndgen.Gaus(P4_min_t[0],sigma);
    P4_min[1] = rndgen.Gaus(P4_min_t[1],sigma);
    P4_min[2] = det2 + B_event::L + B_event::Delta_z;

    GetThetaPhi(P1_plus, P2_plus, theta_plus_in, phi_plus_in);
    GetThetaPhi(P1_min, P2_min, theta_min_in, phi_min_in);
    GetThetaPhi(P3_plus, P4_plus, theta_plus_out, phi_plus_out);
    GetThetaPhi(P3_min, P4_min, theta_min_out, phi_min_out);

    pi_plus_p_reco = 0.299792458*helix::B*B_event::L/(Sqrt(Sin(theta_plus_out)*Sin(theta_plus_out)-Sin(theta_plus_in)*Sin(theta_plus_in)*Cos(phi_plus_in)*Cos(phi_plus_in))*Sign(1.0,helix::B)-Sin(theta_plus_in)*Sin(phi_plus_in));

    pi_min_p_reco = -0.299792458*helix::B*B_event::L/(-Sqrt(Sin(theta_min_out)*Sin(theta_min_out)-Sin(theta_min_in)*Sin(theta_min_in)*Cos(phi_min_in)*Cos(phi_min_in))*Sign(1.0,helix::B)-Sin(theta_min_in)*Sin(phi_min_in));

    // Reconstructed momentum vector
    pi_plus_px_reco = pi_plus_p_reco*Sin(theta_plus_in)*Cos(phi_plus_in);
    pi_plus_py_reco = pi_plus_p_reco*Sin(theta_plus_in)*Sin(phi_plus_in);
    pi_plus_pz_reco = pi_plus_p_reco*Cos(theta_plus_in);

    pi_min_px_reco = pi_min_p_reco*Sin(theta_min_in)*Cos(phi_min_in);
    pi_min_py_reco = pi_min_p_reco*Sin(theta_min_in)*Sin(phi_min_in);
    pi_min_pz_reco = pi_min_p_reco*Cos(theta_min_in);

    p_reco_init[0] = pi_plus_px_reco;
    p_reco_init[1] = pi_plus_py_reco;
    p_reco_init[2] = pi_plus_pz_reco;

    p_reco_init[3] = pi_min_px_reco;
    p_reco_init[4] = pi_min_py_reco;
    p_reco_init[5] = pi_min_pz_reco;

    p_plus_reco_hist->Fill(pi_plus_modp - pi_plus_p_reco);
    p_min_reco_hist->Fill(pi_min_modp - pi_min_p_reco);

    px_plus_reco_hist->Fill(pi_plus_px - pi_plus_px_reco);
    py_plus_reco_hist->Fill(pi_plus_py - pi_plus_py_reco);
    pz_plus_reco_hist->Fill(pi_plus_pz - pi_plus_pz_reco);

    px_min_reco_hist->Fill(pi_min_px - pi_min_px_reco);
    py_min_reco_hist->Fill(pi_min_py - pi_min_py_reco);
    pz_min_reco_hist->Fill(pi_min_pz - pi_min_pz_reco);
    
    K_reco_hist->Fill(K_p - p_reco_init[2] - p_reco_init[5]);

    ////////// ITERATIVE KINEMATIC FIT //////////////

    iter = new IterKinFit;
    iter->Initialize(6,3,p_reco_init, Momentum_Constr_Vector, Momentum_Constr_Derivative, p_reco_err);
    iter->SetThreshold(1E-12);
    
    if(ev_no == 1)
      iter->Minimize(p_reco_final, graph);
    else 
      iter->Minimize(p_reco_final);

    if (ev_no == 1) {
      iter->GetConstraintVector(p_reco_init).Print();
      iter->GetConstraintVector(p_reco_final).Print();      
      iter->PrintResult();
    }

    pi_plus_p_final = Sqrt(Power(p_reco_final[0],2) + Power(p_reco_final[1],2) + Power(p_reco_final[2],2));
    pi_min_p_final = Sqrt(Power(p_reco_final[3],2) + Power(p_reco_final[4],2) + Power(p_reco_final[5],2));

    p_plus_final_hist->Fill(pi_plus_modp - pi_plus_p_final);
    p_min_final_hist->Fill(pi_min_modp - pi_min_p_final);


    px_plus_final_hist->Fill(pi_plus_px - p_reco_final[0]);
    py_plus_final_hist->Fill(pi_plus_py - p_reco_final[1]);
    pz_plus_final_hist->Fill(pi_plus_pz - p_reco_final[2]) ;

    px_min_final_hist->Fill(pi_min_px - p_reco_final[3]);
    py_min_final_hist->Fill(pi_min_py - p_reco_final[4]);
    pz_min_final_hist->Fill(pi_min_pz - p_reco_final[5]);

    K_final_hist->Fill(K_p - p_reco_final[2] - p_reco_final[5]);
    K_final2_hist->Fill(K_p - pi_plus_p_final*Cos(theta_plus_in) - pi_min_p_final*Cos(theta_min_in));
    K_final3_hist->Fill(K_p - Sqrt((Power(p_reco_final[0],2)+Power(p_reco_final[1],2))/Power(Sin(theta_plus_in),2))*Cos(theta_plus_in) - Sqrt((Power(p_reco_final[3],2)+Power(p_reco_final[4],2))/Power(Sin(theta_min_in),2))*Cos(theta_min_in));

    delete iter;
    

    ReadImpr(impr_file, ev_no, K_z, K_p, pi_plus_modp, theta_plus_in_t, phi_plus_in_t, pi_min_modp, theta_min_in_t, phi_min_in_t, P1_plus, P1_min, P2_plus, P2_min);
  }
  while (!impr_file.eof());


  p_plus_reco_hist->Write();
  p_min_reco_hist->Write();
  px_plus_reco_hist->Write();
  py_plus_reco_hist->Write();
  pz_plus_reco_hist->Write();
  px_min_reco_hist->Write();
  py_min_reco_hist->Write();
  pz_min_reco_hist->Write();
  px_plus_final_hist->Write();
  py_plus_final_hist->Write();
  pz_plus_final_hist->Write();
  px_min_final_hist->Write();
  py_min_final_hist->Write();
  pz_min_final_hist->Write();
  p_plus_final_hist->Write();
  p_min_final_hist->Write();
  K_reco_hist->Write();
  K_final_hist->Write();
  K_final2_hist->Write();
  K_final3_hist->Write();
  for (int i = 0; i < 6; i++)
    graph[i].Write();

  root_out->Close();

}

