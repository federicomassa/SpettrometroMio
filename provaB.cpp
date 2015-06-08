#include "B_event.h"
#include <TGraph2D.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace TMath;

void provaB() {

  B_event ev;
  Double_t p1[3] = {0,0,0};
  Double_t p2[3] = {0.004,0.001,1};
  Double_t gamma = 1000;
  Double_t charge = -1;
  
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


  
  cout << fixed << setprecision(8) << "p = " << reco << '\t' << reco2 << endl;
  

}
