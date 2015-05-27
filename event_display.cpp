#include "ReadEvent.h"
#include "ReadDetector.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>

using namespace std;

void event_display(const char* gen_file, const char* det_file, Int_t decay_number) {

  string path = "../Spettrometro_Files/";
  string gen_file_str, det_file_str;
  string str1, str2;

  stringstream ss;
  
  int dec_no; int ev_no_d; double x1_plus; double y1_plus; double z1_plus; double x1_min; double y1_min; double z1_min; double x2_plus; double y2_plus; double z2_plus; double x2_min; double y2_min; double z2_min;

  int ev_no; double K_z; double K_p; double pi_plus_modp; double pi_plus_theta; double pi_plus_phi; double pi_min_modp; double pi_min_theta; double pi_min_phi;

  ss << gen_file;
  ss >> gen_file_str;

  ss.clear();

  ss << det_file;
  ss >> det_file_str;

  ss.clear();

  gen_file_str = path + gen_file_str;
  det_file_str = path + det_file_str;

  ifstream gen_in(gen_file_str.c_str());
  ifstream det_in(det_file_str.c_str());

  for (int i = 0; i < 5; i++) {
    getline(gen_in, str1);
    getline(det_in, str2);
  }

  do {
    ReadEvent(gen_in, ev_no, K_z, K_p, pi_plus_modp, pi_plus_theta, pi_plus_phi, pi_min_modp, pi_min_theta, pi_min_phi);    
  }
  
  while (ev_no != decay_number);

  do {
    ReadDetector(det_in, dec_no, ev_no_d,  x1_plus,  y1_plus,  z1_plus,  x1_min,  y1_min,  z1_min,  x2_plus,  y2_plus,  z2_plus,  x2_min,  y2_min,  z2_min);    
  }
  
  while (ev_no_d != decay_number);

  Double_t x[4];
  Double_t y[4];
  Double_t z[4];

  x[0] = x1_plus;
  x[1] = x2_plus;
  x[2] = x1_min;
  x[3] = x2_min;

  y[0] = y1_plus;
  y[1] = y2_plus;
  y[2] = y1_min;
  y[3] = y2_min;

  z[0] = z1_plus;
  z[1] = z2_plus;
  z[2] = z1_min;
  z[3] = z2_min;

  Double_t decx[1];
  Double_t decy[1];
  Double_t decz[1];

  decx[0] = 0;
  decy[0] = 0;
  decz[0] = K_z;

  Double_t x_t[4];
  Double_t y_t[4];
  Double_t z_t[4];

  x_t[0] = (z1_plus - K_z)*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi);
  x_t[1] = (z2_plus - K_z)*TMath::Tan(pi_plus_theta)*TMath::Cos(pi_plus_phi);
  x_t[2] = (z1_plus - K_z)*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);
  x_t[3] = (z2_plus - K_z)*TMath::Tan(pi_min_theta)*TMath::Cos(pi_min_phi);

  y_t[0] = (z1_plus - K_z)*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
  y_t[1] = (z2_plus - K_z)*TMath::Tan(pi_plus_theta)*TMath::Sin(pi_plus_phi);
  y_t[2] = (z1_plus - K_z)*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);
  y_t[3] = (z2_plus - K_z)*TMath::Tan(pi_min_theta)*TMath::Sin(pi_min_phi);
 
  z_t[0] = z1_plus;
  z_t[1] = z2_plus;
  z_t[2] = z1_min;
  z_t[3] = z2_min;


  TGraph* meas_g_zx = new TGraph(4, z, x);
  // meas_g_zx->GetXaxis()->SetRangeUser(-50,50);
  // meas_g_zx->GetYaxis()->SetRangeUser(-0.2,0.2);
  TGraph* dec_g_zx = new TGraph(1, decz, decx);
  // dec_g_zx->GetXaxis()->SetRangeUser(-50,50);
  // dec_g_zx->GetYaxis()->SetRangeUser(-0.2,0.2);
  TGraph* th_g_zx = new TGraph(4, z_t, x_t);
  // th_g_zx->GetXaxis()->SetRangeUser(-50,50);
  // th_g_zx->GetYaxis()->SetRangeUser(-0.2,0.2);

  meas_g_zx->SetMarkerStyle(kFullCircle);
  meas_g_zx->SetMarkerSize(1);
  meas_g_zx->SetMarkerColor(kRed);

  dec_g_zx->SetMarkerStyle(kFullCircle);
  dec_g_zx->SetMarkerSize(1);
  dec_g_zx->SetMarkerColor(kGreen);

  th_g_zx->SetMarkerStyle(kFullCircle);
  th_g_zx->SetMarkerSize(1);
  th_g_zx->SetMarkerColor(kBlue);

  TCanvas* c1 = new TCanvas();
  c1->SetGrid();
  c1->cd();
  TMultiGraph* mg_zx = new TMultiGraph("mg_zx", "Event Display;z;x");
  mg_zx->Add(meas_g_zx);
  mg_zx->Add(dec_g_zx);
  mg_zx->Add(th_g_zx);
  // mg_zx->GetXaxis()->SetLimits(0,50);
  // mg_zx->GetYaxis()->SetLimits(0,0.2);
  mg_zx->Draw("ap");
  
  ///////////////////////

  TGraph* meas_g_zy = new TGraph(4, z, y);
  // meas_g_zy->GetXaxis()->SetRangeUser(-50,50);
  // meas_g_zy->GetYaxis()->SetRangeUser(-0.2,0.2);
  TGraph* dec_g_zy = new TGraph(1, decz, decy);
  // dec_g_zy->GetXaxis()->SetRangeUser(-50,50);
  // dec_g_zy->GetYaxis()->SetRangeUser(-0.2,0.2);
  TGraph* th_g_zy = new TGraph(4, z_t, y_t);
  // th_g_zy->GetXaxis()->SetRangeUser(-50,50);
  // th_g_zy->GetYaxis()->SetRangeUser(-0.2,0.2);

  meas_g_zy->SetMarkerStyle(kFullCircle);
  meas_g_zy->SetMarkerSize(1);
  meas_g_zy->SetMarkerColor(kRed);

  dec_g_zy->SetMarkerStyle(kFullCircle);
  dec_g_zy->SetMarkerSize(1);
  dec_g_zy->SetMarkerColor(kGreen);

  th_g_zy->SetMarkerStyle(kFullCircle);
  th_g_zy->SetMarkerSize(1);
  th_g_zy->SetMarkerColor(kBlue);

  TCanvas* c2 = new TCanvas();
  c2->SetGrid();
  c2->cd();

  TMultiGraph* mg_zy = new TMultiGraph("mg_zy", "Event Display;z;y");
  mg_zy->Add(meas_g_zy);
  mg_zy->Add(dec_g_zy);
  mg_zy->Add(th_g_zy);
  // mg_zy->GetXaxis()->SetLimits(-50,50);
  // mg_zy->GetYaxis()->SetLimits(-0.2,0.2);
  mg_zy->Draw("ap");


  
  


}
