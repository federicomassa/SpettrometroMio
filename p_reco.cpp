#include "B_event.h"
#include "GetThetaPhi.h"
#include "ReadImpr.h"
#include <fstream>
#include <string>
#include <sstream>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom3.h>

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
  Double_t pi_mass = 0.13957018; // [GeV/c^2]
  B_event ev_plus;
  B_event ev_min;

  int count = 0;
  //HISTOGRAMS
  TH1F* p_plus_reco_hist = new TH1F("p_plus_reco", "Pi+ p Reco; True-Meas", 1000, 1, 0);
  TH1F* p_min_reco_hist = new TH1F("p_min_reco", "Pi- p Reco; True-Meas", 1000, 1, 0);
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
  TFile* root_out = new TFile("p_reco.root","recreate");
  
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

    p_plus_reco_hist->Fill(pi_plus_modp - pi_plus_p_reco);
    p_min_reco_hist->Fill(pi_min_modp - pi_min_p_reco);
    
    ReadImpr(impr_file, ev_no, K_z, K_p, pi_plus_modp, theta_plus_in_t, phi_plus_in_t, pi_min_modp, theta_min_in_t, phi_min_in_t, P1_plus, P1_min, P2_plus, P2_min);
  }
  while (!impr_file.eof());


  p_plus_reco_hist->Write();
  p_min_reco_hist->Write();
  root_out->Close();

}

