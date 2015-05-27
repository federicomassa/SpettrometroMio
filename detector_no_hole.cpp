#include "ReadEvent.h"
#include <fstream>
#include <string>
#include <TMath.h>
#include <TRandom3.h>
#include <iomanip>
#include <ctype.h>

void detector(const char* in_filename = "events.dat", const char* out_filename = "measures.dat") {
  
  string path = "../Spettrometro_Files/";
  string in_filename_str, out_filename_str;
  stringstream ss;
  ss << in_filename;
  ss >> in_filename_str;
  ss.clear();
  
  ss << out_filename;
  ss >> out_filename_str;
  
  string complete_in_filename = path + in_filename_str;
  string complete_out_filename = path + out_filename_str;
  
  ifstream in(complete_in_filename.c_str()); //input file (events)
  ofstream detector_out(complete_out_filename.c_str()); //output file (detector response)
  string str;
  string ev_str;
  getline(in,str);
  detector_out << str << endl; //Time info
  getline(in,str);
  detector_out << str << endl; //N. of events 
  ev_str = str.substr(10);
  int imax = atoi(ev_str.c_str()); //Save n. of events
  detector_out << '\n';
  getline(in,str);
  getline(in,str);
  getline(in,str);
  
  ///// EVENT VARIABLES
  int ev_no = 0, dec_no;
  double K_z, K_p;
  double pi_plus_modp, pi_plus_theta, pi_plus_phi;
  double pi_min_modp, pi_min_theta, pi_min_phi;
  

  ///// DETECTOR VARIABLES
  TRandom3 rndgen; // Pseudo-Random numbers generator
  double P1_plus_t[3], P1_min_t[3]; //First detector true points
  double P1_plus_m[3], P1_min_m[3]; //First detector measured points
  double P2_plus_t[3], P2_min_t[3]; //Second detector true points
  double P2_plus_m[3], P2_min_m[3]; //Second detector measured points
  double z1 = 25; //First detector z position [m]
  double z2 = 35; //Second detector z position [m]
  P1_plus_t[2] = z1;
  P1_min_t[2] = z1;
  P2_plus_t[2] = z2;
  P2_min_t[2] = z2;
  double sigma_x = 0.001, sigma_y = 0.001;  //detector resolution [m]
  double d1, d2; //event distance from detectors
  
  ///// Write legend /////////////////////
  detector_out << "Dec_no" << '\t' << "Ev_no" << '\t' << "x1+" << '\t' << "y1+" << '\t' << "z1+" << '\t' << "x1-" << '\t' << "y1-" << '\t' << "z1-" << '\t' << "x2+" << '\t' << "y2+" << '\t' << "z2+" << '\t' << "x2-" << '\t' << "y2-" << '\t' << "z2-" << '\n' << '\n';

  // int k = 1;

  ////////// EVENT CYCLE ////////////
  ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);


  while (!in.eof()) {
    // k++;
    if (dec_no % (int(double(imax)/20)) == 0) cout << double(dec_no)/double(imax)*100 << "% completed..." << endl;

    // cout << ev_no << '\t' << K_z << '\t' << K_p << '\t' << pi_plus_modp << '\t' << pi_plus_theta << '\t' << pi_plus_phi << '\t' << pi_min_modp << '\t' << pi_min_theta << '\t' << pi_min_phi << endl;
    d1 = z1 - K_z;
    d2 = z2 - K_z;
    if (d1 >= 0) {
      ev_no++;
      P1_plus_t[0] = d1*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi); //x1,y1 calculation given z1, theta, phi
      P1_min_t[0] = d1*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);
      P1_plus_t[1] = d1*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
      P1_min_t[1] = d1*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);
      P2_plus_t[0] = d2*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi); //x2,y2 calculation given z2, theta, phi
      P2_min_t[0] = d2*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);
      P2_plus_t[1] = d2*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
      P2_min_t[1] = d2*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);

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

      detector_out << fixed << setprecision(0) << dec_no << '\t'  << ev_no << '\t' << setprecision(5) << P1_plus_m[0] << '\t' << P1_plus_m[1] << '\t' << P1_plus_m[2] << '\t' << P1_min_m[0] << '\t' << P1_min_m[1] << '\t' << P1_min_m[2] << '\t' << P2_plus_m[0] << '\t' << P2_plus_m[1] << '\t' << P2_plus_m[2] << '\t' << P2_min_m[0] << '\t' << P2_min_m[1] << '\t' << P2_min_m[2] << endl;
      
    }
    ReadEvent(in,dec_no,K_z,K_p,pi_plus_modp,pi_plus_theta,pi_plus_phi,pi_min_modp,pi_min_theta,pi_min_phi);
  }
  in.close();
  detector_out.close();
}
