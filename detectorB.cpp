#include "ReadEvent.h"
#include "B_event.h"
#include <fstream>
#include <string>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom3.h>
#include <iomanip>
#include <ctype.h>

void detectorB(const char* in_filename = "events", const char* out_filename = "measuresB") {
  
  //// FILES WRITING

  string path = "../Spettrometro_Files/";
  string dat_ext = ".dat";
  string root_ext = ".root";
  string in_filename_str, out_filename_str;
  stringstream ss;
  ss << in_filename;
  ss >> in_filename_str;
  ss.clear();
  
  ss << out_filename;
  ss >> out_filename_str;
  ss.clear();

  string complete_in_filename = path + in_filename_str + dat_ext;
  string complete_out_filename = path + out_filename_str + dat_ext;
  string complete_root_filename = path + out_filename_str + root_ext;

  ifstream in(complete_in_filename.c_str()); //input file (events)
  ofstream detector_out(complete_out_filename.c_str()); //output file (detector response)
  TFile* root_out = new TFile(complete_root_filename.c_str(), "RECREATE");

  //// GENERAL EVENT INFO
  Double_t pi_mass = 0.13957018;
  string str;
  string ev_str;
  getline(in,str);
  detector_out << str << endl; //Time info
  getline(in,str);
  detector_out << str << endl; //N. of events 
  ev_str = str.substr(10);
  Int_t imax = atoi(ev_str.c_str()); //Save n. of events
  detector_out << '\n';
  getline(in,str);
  getline(in,str);
  getline(in,str);
  
  ///// EVENT VARIABLES
  Int_t ev_no = 0, dec_no;
  Double_t K_z, K_p;
  Double_t pi_plus_modp, pi_plus_theta, pi_plus_phi;
  Double_t pi_min_modp, pi_min_theta, pi_min_phi;
  B_event mag_event_plus, mag_event_min;

  ///// DETECTOR VARIABLES
  TRandom3 rndgen; // Pseudo-Random numbers generator

  Double_t P1_plus_t[3], P1_min_t[3]; //First detector true points
  Double_t P1_plus_m[3], P1_min_m[3]; //First detector measured points
  Double_t P2_plus_t[3], P2_min_t[3]; //Second detector true points
  Double_t P2_plus_m[3], P2_min_m[3]; //Second detector measured points

  Double_t P3_plus_t[3], P3_min_t[3]; //Third detector true points
  Double_t P3_plus_m[3], P3_min_m[3]; //Third detector measured points
  Double_t P4_plus_t[3], P4_min_t[3]; //Fourth detector true points
  Double_t P4_plus_m[3], P4_min_m[3]; //Fourth detector measured points



  Double_t z1 = 50; //First detector z position [m]
  Double_t z2 = 60; //Second detector z position [m]
  
  // Last two detector coordinates written in class B_event
  // Double_t z3 = z2 + B_event::L; //Third detector z position [m]
  // Double_t z4 = z3 + B_event::Delta_z; //Fourth detector z position [m]

  Double_t R_int1 = 0.084; //internal and external detector radius [m] (only R_int1 is important)
  // Double_t R_ext1 = 0.40; //Check values if needed. Theta_max = 0.007485, Max(theta_minimum) = 0.004195
  // Double_t R_int2 = 0.20;
  // Double_t R_ext2 = 0.60;
  P1_plus_t[2] = z1;
  P1_min_t[2] = z1;
  P2_plus_t[2] = z2;
  P2_min_t[2] = z2;
  Double_t sigma_x = 0.001, sigma_y = 0.001;  //detector resolution [m]
  Double_t d1, d2; //event distance from detectors

 

  ////// HISTOGRAMS
  TH1F* pi_plus_theta_hist = new TH1F("pi_plus_theta_hist", "Pi+ Theta Accepted Hist;Theta;#", 1000,1,0);
  TH1F* pi_min_theta_hist = new TH1F("pi_min_theta_hist", "Pi- Theta Accepted Hist;Theta;#", 1000,1,0);
  TH1F* pi_plus_phi_hist = new TH1F("pi_plus_phi_hist", "Pi+ Phi Accepted Hist;Phi;#", 1000,1,0);
  TH1F* pi_min_phi_hist = new TH1F("pi_min_phi_hist", "Pi- Phi Accepted Hist;Phi;#", 1000,1,0);
  
  ///// Write legend /////////////////////
  detector_out << "Dec_no" << '\t' << "Ev_no" << '\t' << "x1+" << '\t' << "y1+" << '\t' << "z1+" << '\t' << "x1-" << '\t' << "y1-" << '\t' << "z1-" << '\t' << "x2+" << '\t' << "y2+" << '\t' << "z2+" << '\t' << "x2-" << '\t' << "y2-" << '\t' << "z2-" << '\t' << "x3+" << '\t' << "y3+" << '\t' << "z3+" << '\t' << "x3-" << '\t' << "y3-" << '\t' << "z3-" << '\t' << "x4+" << '\t' << "y4+" << '\t' << "z4+" << '\t' << "x4-" << '\t' << "y4-" << '\t' << "z4-" << '\n' << '\n';

  // int k = 1;

  ////////// EVENT CYCLE ////////////
  ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);

  while (!in.eof()) {
    // k++;
    if (dec_no % (int(Double_t(imax)/100*5)) == 0) cout << int(Double_t(dec_no)/Double_t(imax)*100) << "% completed..." << endl;

    // cout << ev_no << '\t' << K_z << '\t' << K_p << '\t' << pi_plus_modp << '\t' << pi_plus_theta << '\t' << pi_plus_phi << '\t' << pi_min_modp << '\t' << pi_min_theta << '\t' << pi_min_phi << endl;
    d1 = z1 - K_z;
    d2 = z2 - K_z;
    if (d1 >= 0) {
      
      P1_plus_t[0] = d1*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi); //x1,y1 calculation given z1, theta, phi
      P1_min_t[0] = d1*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);
      P1_plus_t[1] = d1*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
      P1_min_t[1] = d1*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);
      P2_plus_t[0] = d2*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi); //x2,y2 calculation given z2, theta, phi
      P2_min_t[0] = d2*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);
      P2_plus_t[1] = d2*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
      P2_min_t[1] = d2*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);

      //Apply magnetic field
      mag_event_plus.SetB_event(P1_plus_t, P2_plus_t, 1.0, TMath::Sqrt(1+TMath::Power(pi_plus_modp/pi_mass,2)));
      mag_event_min.SetB_event(P1_min_t, P2_min_t, 1.0, TMath::Sqrt(1+TMath::Power(pi_min_modp/pi_mass,2)));

      mag_event_plus.GetP3(P3_plus_t);
      mag_event_min.GetP3(P3_min_t);
      
      mag_event_plus.GetP4(P4_plus_t);
      mag_event_min.GetP4(P4_min_t);

      // Acceptance
      if (P1_plus_t[0]*P1_plus_t[0] + P1_plus_t[1]*P1_plus_t[1] < R_int1*R_int1 ||
	  P1_min_t[0]*P1_min_t[0] + P1_min_t[1]*P1_min_t[1] < R_int1*R_int1
	  ) {
	ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);
	continue;}

      ev_no++; //Event has been accepted

      pi_plus_theta_hist->Fill(pi_plus_theta);
      pi_min_theta_hist->Fill(pi_min_theta);
      pi_plus_phi_hist->Fill(pi_plus_phi);
      pi_min_phi_hist->Fill(pi_min_phi);

      // Detector resolution
      //     Before magnetic field
      P1_plus_m[0] = rndgen.Gaus(P1_plus_t[0],sigma_x);
      P1_plus_m[1] = rndgen.Gaus(P1_plus_t[1],sigma_y);
      P1_plus_m[2] = P1_plus_t[2];

      P1_min_m[0] = rndgen.Gaus(P1_min_t[0],sigma_x);
      P1_min_m[1] = rndgen.Gaus(P1_min_t[1],sigma_y);
      P1_min_m[2] = P1_min_t[2];

      P2_plus_m[0] = rndgen.Gaus(P2_plus_t[0],sigma_x);
      P2_plus_m[1] = rndgen.Gaus(P2_plus_t[1],sigma_y);
      P2_plus_m[2] = P2_plus_t[2];

      P2_min_m[0] = rndgen.Gaus(P2_min_t[0],sigma_x);
      P2_min_m[1] = rndgen.Gaus(P2_min_t[1],sigma_y);
      P2_min_m[2] = P2_min_t[2];

      //     After magnetic field
      P3_plus_m[0] = rndgen.Gaus(P3_plus_t[0],sigma_x);
      P3_plus_m[1] = rndgen.Gaus(P3_plus_t[1],sigma_y);
      P3_plus_m[2] = P3_plus_t[2];

      P3_min_m[0] = rndgen.Gaus(P3_min_t[0],sigma_x);
      P3_min_m[1] = rndgen.Gaus(P3_min_t[1],sigma_y);
      P3_min_m[2] = P3_min_t[2];

      P4_plus_m[0] = rndgen.Gaus(P4_plus_t[0],sigma_x);
      P4_plus_m[1] = rndgen.Gaus(P4_plus_t[1],sigma_y);
      P4_plus_m[2] = P4_plus_t[2];

      P4_min_m[0] = rndgen.Gaus(P4_min_t[0],sigma_x);
      P4_min_m[1] = rndgen.Gaus(P4_min_t[1],sigma_y);
      P4_min_m[2] = P4_min_t[2];

      

      //if measure falls outside the detector (even if the theoretical one
      //was inside, bring the measure to the closest point of the detector
      if (P1_plus_m[0]*P1_plus_m[0] + P1_plus_m[1]*P1_plus_m[1] < R_int1*R_int1)
	{
	  P1_plus_m[0] = R_int1*P1_plus_m[0]/
	    TMath::Sqrt(P1_plus_m[0]*P1_plus_m[0] + P1_plus_m[1]*P1_plus_m[1]);

	  P1_plus_m[1] = R_int1*P1_plus_m[1]/
	    TMath::Sqrt(P1_plus_m[0]*P1_plus_m[0] + P1_plus_m[1]*P1_plus_m[1]);

	   
	}

      if (P1_min_m[0]*P1_min_m[0] + P1_min_m[1]*P1_min_m[1] < R_int1*R_int1)
	{
	  P1_min_m[0] = R_int1*P1_min_m[0]/
	    TMath::Sqrt(P1_min_m[0]*P1_min_m[0] + P1_min_m[1]*P1_min_m[1]);

	  P1_min_m[1] = R_int1*P1_min_m[1]/
	    TMath::Sqrt(P1_min_m[0]*P1_min_m[0] + P1_min_m[1]*P1_min_m[1]);
	}

      // Write detector response
      detector_out << fixed << setprecision(0) << dec_no << '\t'  << ev_no << '\t' << setprecision(5) << P1_plus_m[0] << '\t' << P1_plus_m[1] << '\t' << P1_plus_m[2] << '\t' << P1_min_m[0] << '\t' << P1_min_m[1] << '\t' << P1_min_m[2] << '\t' << P2_plus_m[0] << '\t' << P2_plus_m[1] << '\t' << P2_plus_m[2] << '\t' << P2_min_m[0] << '\t' << P2_min_m[1] << '\t' << P2_min_m[2] << '\t' << P3_plus_m[0] << '\t' << P3_plus_m[1] << '\t' << P3_plus_m[2] << '\t' << P3_min_m[0] << '\t' << P3_min_m[1] << '\t' << P3_min_m[2] << '\t' << P3_plus_m[0] << '\t' << P4_plus_m[1] << '\t' << P4_plus_m[2] << '\t' << P4_min_m[0] << '\t' << P4_min_m[1] << '\t' << P4_min_m[2] << endl;
      
    }
    ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);
  }

  pi_plus_theta_hist->Write();
  pi_min_theta_hist->Write();
  pi_plus_phi_hist->Write();
  pi_min_phi_hist->Write();

  in.close();
  detector_out.close();
  root_out->Close();
}
