#include "ReadEvent.h"
#include "B_event_approx.h"
#include <fstream>
#include <string>
#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom3.h>
#include <iomanip>
#include <ctype.h>
#include <TNtupleD.h>

void detector(const char* in_filename = "events", const char* out_ntuple = "measures") {
  
  //// FILES WRITING

  string path = "../Spettrometro_Files/";
  string dat_ext = ".dat";
  string root_ext = ".root";
  string in_filename_str, out_ntuple_str;
  stringstream ss;
  ss << in_filename;
  ss >> in_filename_str;
  ss.clear();

  ss << out_ntuple;
  ss >> out_ntuple_str;

  string complete_in_filename = path + in_filename_str + dat_ext;
  string complete_root_filename = path + out_ntuple_str + root_ext;

  string ntuple_title = "Generated event and detector response";
  
  string branch_list = "K_z:K_p:p_plus:p_min:";
  branch_list += "x1_plus_t:";
  branch_list += "y1_plus_t:";
  branch_list += "x2_plus_t:";
  branch_list += "y2_plus_t:";
  branch_list += "x3_plus_t:";
  branch_list += "y3_plus_t:";
  branch_list += "x4_plus_t:";
  branch_list += "y4_plus_t:";
  branch_list += "x1_min_t:";
  branch_list += "y1_min_t:";
  branch_list += "x2_min_t:";
  branch_list += "y2_min_t:";
  branch_list += "x3_min_t:";
  branch_list += "y3_min_t:";
  branch_list += "x4_min_t:";
  branch_list += "y4_min_t:";
  branch_list += "x1_plus_m:";
  branch_list += "y1_plus_m:";
  branch_list += "x2_plus_m:";
  branch_list += "y2_plus_m:";
  branch_list += "x3_plus_m:";
  branch_list += "y3_plus_m:";
  branch_list += "x4_plus_m:";
  branch_list += "y4_plus_m:";
  branch_list += "x1_min_m:";
  branch_list += "y1_min_m:";
  branch_list += "x2_min_m:";
  branch_list += "y2_min_m:";
  branch_list += "x3_min_m:";
  branch_list += "y3_min_m:";
  branch_list += "x4_min_m:";
  branch_list += "y4_min_m";


  ifstream in(complete_in_filename.c_str()); //input file (events)
  TNtupleD nt_det(out_ntuple, ntuple_title.c_str(), branch_list.c_str());
  TFile* root_out = new TFile(complete_root_filename.c_str(), "RECREATE");

  //// GENERAL EVENT INFO
  string str;
  string ev_str;
  getline(in,str);
  getline(in,str);
  ev_str = str.substr(10);
  int imax = atoi(ev_str.c_str()); //Save n. of events
  getline(in,str);
  getline(in,str);
  getline(in,str);
  
  ///// EVENT VARIABLES
  int ev_no = 0, dec_no;
  Double_t K_z, K_p;
  Double_t pi_plus_modp, pi_plus_theta, pi_plus_phi;
  Double_t pi_min_modp, pi_min_theta, pi_min_phi;
  B_event_approx *ev_plus, *ev_min;

  Double_t* nt_fill = new Double_t[36];

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
  Double_t R_int1 = 0.084; //internal and external detector radius [m] (only R_int1 is important)
  // double R_ext1 = 0.40; //Check values if needed. Theta_max = 0.007485, Max(theta_minimum) = 0.004195
  // double R_int2 = 0.20;
  // double R_ext2 = 0.60;
  P1_plus_t[2] = z1;
  P1_min_t[2] = z1;
  P2_plus_t[2] = z2;
  P2_min_t[2] = z2;
  Double_t sigma_x = 0.001, sigma_y = 0.001;  //detector resolution [m]
  Double_t d1, d2; //event distance from detectors

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

      

      // Acceptance
      if (P1_plus_t[0]*P1_plus_t[0] + P1_plus_t[1]*P1_plus_t[1] < R_int1*R_int1 ||
	  P1_min_t[0]*P1_min_t[0] + P1_min_t[1]*P1_min_t[1] < R_int1*R_int1
	  ) {
	ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi); 
	continue;}

      ev_no++; //Event has been accepted

      ev_plus = new B_event_approx;
      ev_min = new B_event_approx;
      
      ev_plus->SetB_event_approx(P1_plus_t, P2_plus_t, 1.0, pi_plus_modp);
      ev_min->SetB_event_approx(P1_min_t, P2_min_t, -1.0, pi_min_modp);

      ev_plus->GetP3(P3_plus_t);
      ev_plus->GetP4(P4_plus_t);

      ev_min->GetP3(P3_min_t);
      ev_min->GetP4(P4_min_t);

      // Detector resolution
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

      nt_fill[0] = K_z;
      nt_fill[1] = K_p;
      nt_fill[2] = pi_plus_modp;
      nt_fill[3] = pi_min_modp;
      nt_fill[4] = P1_plus_t[0];
      nt_fill[5] = P1_plus_t[1];
      nt_fill[6] = P2_plus_t[0];
      nt_fill[7] = P2_plus_t[1];
      nt_fill[8] = P3_plus_t[0];
      nt_fill[9] = P3_plus_t[1];
      nt_fill[10] = P4_plus_t[0];
      nt_fill[11] = P4_plus_t[1];
      nt_fill[12] = P1_min_t[0];
      nt_fill[13] = P1_min_t[1];
      nt_fill[14] = P2_min_t[0];
      nt_fill[15] = P2_min_t[1];
      nt_fill[16] = P3_min_t[0];
      nt_fill[17] = P3_min_t[1];
      nt_fill[18] = P4_min_t[0];
      nt_fill[19] = P4_min_t[1];
      nt_fill[20] = P1_plus_m[0];
      nt_fill[21] = P1_plus_m[1];
      nt_fill[22] = P2_plus_m[0];
      nt_fill[23] = P2_plus_m[1];
      nt_fill[24] = P3_plus_m[0];
      nt_fill[25] = P3_plus_m[1];
      nt_fill[26] = P4_plus_m[0];
      nt_fill[27] = P4_plus_m[1];
      nt_fill[28] = P1_min_m[0];
      nt_fill[29] = P1_min_m[1];
      nt_fill[30] = P2_min_m[0];
      nt_fill[31] = P2_min_m[1];
      nt_fill[32] = P3_min_m[0];
      nt_fill[33] = P3_min_m[1];
      nt_fill[34] = P4_min_m[0];
      nt_fill[35] = P4_min_m[1];

      delete ev_plus;
      delete ev_min;

  
    }
    ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);
    nt_det.Fill(nt_fill);

    

  }

  nt_det.Write();
  in.close();
  root_out->Close();
}
