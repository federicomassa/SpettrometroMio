#include "ReadEvent.h"
#include "B_event.h"
#include <fstream>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace TMath;

void error_omegaT() {

  ifstream in_file("../Spettrometro_Files/events.dat");
  TFile* root_out = new TFile("error_omegaT.root", "recreate");
  string str;

  TH1F* diff_plus = new TH1F("diff_plus", "Pi+ difference in OmegaT; True-Approx", 100, 1,0);
  TH1F* diff_min = new TH1F("diff_min", "Pi- difference in OmegaT; True-Approx", 100, 1,0);
  TH1F* cos_arg_plus = new TH1F("cos_arg_plus", "omegaT + atan(tan(theta)sin(phi))", 100, 1,0);
  TH1F* cos_arg_min = new TH1F("cos_arg_min", "omegaT + atan(tan(theta)sin(phi))", 100, 1,0);
  TH1F* pi_plus_reco_hist = new TH1F("pi_plus_reco_hist", "Pi+ reco", 10000, 1,0);
  TH1F* pi_min_reco_hist = new TH1F("pi_min_reco_hist", "Pi- reco", 10000, 1,0);
  TH1F* z_dec_hist = new TH1F("z_dec_hist", "z dec Hist", 1000, 1, 0);
  TH1F* K_p_hist = new TH1F("K_p_hist", "z dec Hist", 1000, 1, 0);
  TH1F* pi_plus_modp_hist = new TH1F("pi_plus_modp_hist", "z dec Hist", 1000, 1, 0);
  TH1F* pi_plus_theta_hist = new TH1F("pi_plus_theta_hist", "z dec Hist", 1000, 1, 0);
  TH1F* pi_plus_phi_hist = new TH1F("pi_plus_phi_hist", "z dec Hist", 1000, 1, 0);


  getline(in_file,str);
  getline(in_file,str);
  getline(in_file,str);
  getline(in_file,str);
  getline(in_file,str);

  Int_t ev_no;
  Double_t z_dec, K_p;
  Double_t pi_plus_modp, pi_plus_theta, pi_plus_phi;
  Double_t pi_min_modp, pi_min_theta, pi_min_phi;
  B_event m_event_plus, m_event_min;

  Double_t P1_plus[3], P1_min[3],
    P2_plus[3], P2_min[3],
    P3_plus[3], P3_min[3],
    P4_plus[3], P4_min[3];
  
  Double_t gamma_plus, beta_plus; 
  Double_t gamma_min, beta_min;
  Double_t init_speed_plus[3], init_speed_min[3];
  Double_t pi_mass = 0.13957018;
  Double_t omega_t_plus, radius_plus;
  Double_t omega_t_min, radius_min;
  Double_t theta_plus_out, phi_plus_out, theta_min_out, phi_min_out;
  Double_t pi_plus_reco, pi_min_reco;
  Double_t det1, det2, det3, det4;
  Double_t R_int1 = 0.084;
  
  int count = 0;
  //Coordinates of detectors
  det1 = 50;
  det2 = 60;
  det3 = det2 + B_event::L;
  det4 = det3 + B_event::Delta_z;

  ReadEvent(in_file, ev_no, z_dec, K_p, pi_plus_modp, pi_plus_theta, pi_plus_phi,pi_min_modp, pi_min_theta, pi_min_phi);
 
  do {

 

    z_dec_hist->Fill(z_dec);
    K_p_hist->Fill(K_p);
    pi_plus_modp_hist->Fill(pi_plus_modp);
    pi_plus_theta_hist->Fill(pi_plus_theta);
    pi_plus_phi_hist->Fill(pi_plus_phi);

    if (z_dec > det1) {ReadEvent(in_file, ev_no, z_dec, K_p, pi_plus_modp, pi_plus_theta, pi_plus_phi,pi_min_modp, pi_min_theta, pi_min_phi); continue;}

    //Set points 
    P1_plus[0] = (det1-z_dec)*Tan(pi_plus_theta)*Cos(pi_plus_phi);
    P1_plus[1] = (det1-z_dec)*Tan(pi_plus_theta)*Sin(pi_plus_phi);
    P1_plus[2] = det1;

    P1_min[0] = (det1-z_dec)*Tan(pi_min_theta)*Cos(pi_min_phi);
    P1_min[1] = (det1-z_dec)*Tan(pi_min_theta)*Sin(pi_min_phi);
    P1_min[2] = det1;

    P2_plus[0] = (det2-z_dec)*Tan(pi_plus_theta)*Cos(pi_plus_phi);
    P2_plus[1] = (det2-z_dec)*Tan(pi_plus_theta)*Sin(pi_plus_phi);
    P2_plus[2] = det2;

    P2_min[0] = (det2-z_dec)*Tan(pi_min_theta)*Cos(pi_min_phi);
    P2_min[1] = (det2-z_dec)*Tan(pi_min_theta)*Sin(pi_min_phi);
    P2_min[2] = det2;
    ///////////////////////////

    if (P1_plus[0]*P1_plus[0] + P1_plus[1]*P1_plus[1] < R_int1*R_int1 ||
	  P1_min[0]*P1_min[0] + P1_min[1]*P1_min[1] < R_int1*R_int1
	  ) {
	ReadEvent(in_file,ev_no,z_dec,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi); 
	continue;}
 
    count += 1;

    // cout << ev_no << '\t' << z_dec << '\t' << K_p << '\t' << pi_plus_modp << '\t' << pi_plus_theta << '\t' << pi_plus_phi << '\t' << pi_min_modp << '\t' << pi_min_theta << '\t' << pi_min_phi << endl;

    gamma_plus = Sqrt(Power(pi_plus_modp/pi_mass,2)+1);
    gamma_min = Sqrt(Power(pi_min_modp/pi_mass,2)+1);

    // cout << gamma_plus << '\t' << gamma_min << endl;

    beta_plus = Sqrt(gamma_plus*gamma_plus-1)/gamma_plus;
    beta_min = Sqrt(gamma_min*gamma_min-1)/gamma_min;

    // cout << beta_plus << '\t' << beta_min << endl;

    init_speed_plus[0] = beta_plus*Sin(pi_plus_theta)*Cos(pi_plus_phi);
    init_speed_plus[1] = beta_plus*Sin(pi_plus_theta)*Sin(pi_plus_phi);
    init_speed_plus[2] = beta_plus*Cos(pi_plus_theta);

    init_speed_min[0] = beta_min*Sin(pi_min_theta)*Cos(pi_min_phi);
    init_speed_min[1] = beta_min*Sin(pi_min_theta)*Sin(pi_min_phi);
    init_speed_min[2] = beta_min*Cos(pi_min_theta);


    m_event_plus.SetB_event(P1_plus, P2_plus, 1.0, gamma_plus);
    m_event_min.SetB_event(P1_min, P2_min, -1.0, gamma_min);

    theta_plus_out = m_event_plus.GetTheta_out();
    theta_min_out = m_event_min.GetTheta_out();
    phi_plus_out = m_event_plus.GetPhi_out();
    phi_min_out = m_event_min.GetPhi_out();

    // cout << theta_plus_out << '\t' << phi_plus_out << endl;
    // cout << theta_min_out << '\t' << phi_min_out << endl;

    omega_t_plus = m_event_plus.GetOmegaT();
    omega_t_min = m_event_min.GetOmegaT();

    pi_plus_reco = 0.299792458*1.00*B_event::L/(TMath::Sqrt(TMath::Sin(theta_plus_out)*TMath::Sin(theta_plus_out)-TMath::Sin(pi_plus_theta)*TMath::Sin(pi_plus_theta)*TMath::Cos(pi_plus_phi)*TMath::Cos(pi_plus_phi))*TMath::Sign(1.0,1.00)-TMath::Sin(pi_plus_theta)*TMath::Sin(pi_plus_phi));

    pi_min_reco = -0.299792458*1.00*B_event::L/(-Sqrt(Sin(theta_min_out)*Sin(theta_min_out)-Sin(pi_min_theta)*Sin(pi_min_theta)*Cos(pi_min_phi)*Cos(pi_min_phi))*Sign(1.0,1.00)-Sin(pi_min_theta)*Sin(pi_min_phi));

    // cout << pi_plus_reco << '\t' << pi_min_reco << endl;

    diff_plus->Fill((omega_t_plus-ASin(B_event::L/m_event_plus.GetRadius()))/omega_t_plus);
    diff_min->Fill((omega_t_min+ASin(B_event::L/m_event_min.GetRadius()))/omega_t_min);
    cos_arg_plus->Fill(omega_t_plus + ATan(Tan(pi_plus_theta)*Sin(pi_plus_phi)));
    cos_arg_min->Fill(omega_t_min + ATan(Tan(pi_min_theta)*Sin(pi_min_phi)));
    
    pi_plus_reco_hist->Fill(pi_plus_reco-pi_plus_modp);
    pi_min_reco_hist->Fill(pi_min_reco-pi_min_modp);

    ReadEvent(in_file, ev_no, z_dec, K_p, pi_plus_modp, pi_plus_theta, pi_plus_phi,pi_min_modp, pi_min_theta, pi_min_phi);

  }

  while(!in_file.eof());

  diff_plus->Write();
  diff_min->Write();
  cos_arg_plus->Write();
  cos_arg_min->Write();
  pi_plus_reco_hist->Write();
  pi_min_reco_hist->Write();
  z_dec_hist->Write();
  K_p_hist->Write();
  pi_plus_modp_hist->Write();
  pi_plus_theta_hist->Write();
  pi_plus_phi_hist->Write();

  root_out->Close();
}
