#include "B_event.h"
#include <TGraph2D.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>

using namespace std;
using namespace TMath;

void provaB() {

  cout <<  fixed << setprecision(10);

  B_event ev;
  Double_t p1[3] = {0,0,0};
  Double_t p2[3] = {0.1,0.007,10};
  Double_t p3[3];
  Double_t p4[3];
  Double_t gamma = 70;
  Double_t charge = 1;
  
  TGraph2D* g;
  Double_t theta_in, theta_out, phi_in, phi_out;

  ev.SetB_event(p1, p2, charge, gamma);
  g = ev.GetEventGraph();
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(7);
  g->Draw("P");

  theta_in = ev.GetTheta_in();
  theta_out = ev.GetTheta_out();
  phi_in = ev.GetPhi_in();
  phi_out = ev.GetPhi_out();
  
  cout << "OMEGA_T" << ev.GetTimeInterval()*ev.GetOmega() << endl;
  cout << "TI: " << ev.GetTimeInterval() << endl;
  cout << "RAD: " << ev.GetRadius() << endl;

  ev.GetP3(p3);
  ev.GetP4(p4);

  Double_t x3_t = p3[0];
  Double_t y3_t = p3[1];
  Double_t x4_t = p4[0];
  Double_t y4_t = p4[1];

  Double_t x3_p = p1[0] + TMath::Tan(theta_in)*Cos(phi_in)*(p2[2]+B_event::L);
  Double_t x4_p = p1[0] + TMath::Tan(theta_in)*Cos(phi_in)*(p2[2]+B_event::L+B_event::Delta_z);

  cout << "x3_t: " << p3[0] << endl;
  cout << "x3_p: " << x3_p << endl;

  cout << "x4_t: " << p4[0] << endl;
  cout << "x4_p: " << x4_p << endl;
  
  Double_t reco, reco2;
  Double_t b;
  b = ev.GetTheta_in()*TMath::Sin(ev.GetPhi_in())/TMath::Cos(ev.GetTheta_in());
  Double_t a;
  a = 0.5;
  Double_t c;
  c = (TMath::Cos(ev.GetTheta_out()) - TMath::Cos(ev.GetTheta_in()))/TMath::Cos(ev.GetTheta_in());

  reco = (-b + TMath::Sqrt(b*b-4*a*c))/(2*a);
  cout << "RECO: " << reco << endl;

  reco = TMath::Sqrt(TMath::Power(TMath::Sin(theta_out),2)-TMath::Power(TMath::Sin(theta_in)*TMath::Cos(phi_in),2))*TMath::Sign(1.0,charge)-TMath::Sin(theta_in)*TMath::Sin(phi_in);

  reco = charge*1*B_event::L/reco;
  reco = reco*0.299792458; //UM

  reco2 = 0.299792458*charge*helix::B*B_event::L/(Sqrt(Sin(theta_out)*Sin(theta_out)-Sin(theta_in)*Sin(theta_in)*Cos(phi_in)*Cos(phi_in))*Sign(1.0,charge*helix::B)-Sin(theta_in)*Sin(phi_in));


  
  cout << "p = " << reco << '\t' << reco2 << endl;
  

}
