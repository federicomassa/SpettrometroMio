#include "B_event.h"
#include "B_event_approx.h"
#include <TGraph2D.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <TCanvas.h>

using namespace std;
using namespace TMath;

void provaB() {

  cout <<  fixed << setprecision(10);

  B_event ev;
  B_event_approx ev_approx;
  Double_t p1[3] = {0,0,0};
  Double_t p2[3] = {0.08,0.02,10};
  Double_t p3[3];
  Double_t p4[3];
  Double_t p3_approx[3];
  Double_t p4_approx[3];
  Double_t gamma = 70;
  Double_t charge = -1;
  Double_t pi_mass = 0.13957018;
  Double_t p = pi_mass*gamma*TMath::Sqrt(gamma*gamma-1)/gamma;

  TGraph2D* g;
  TGraph2D* g_approx;
  Double_t theta_in, theta_out, phi_in, phi_out;

  ev.SetB_event(p1, p2, charge, gamma);
  ev_approx.SetB_event_approx(p1, p2, charge, p);
  g = ev.GetEventGraph();
  g_approx = ev_approx.GetEventGraph();

  TCanvas* c1 = new TCanvas();
  c1->cd();
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(7);
  g->Draw("P");

  TCanvas* c2 = new TCanvas();
  c2->cd();
  g_approx->SetMarkerColor(kRed);
  g_approx->SetMarkerStyle(7);
  g_approx->Draw("P");


  theta_in = ev.GetTheta_in();
  theta_out = ev.GetTheta_out();
  phi_in = ev.GetPhi_in();
  phi_out = ev.GetPhi_out();
  
  cout << "OMEGA_T" << ev.GetTimeInterval()*ev.GetOmega() << endl;
  cout << "TI: " << ev.GetTimeInterval() << endl;
  cout << "RAD: " << ev.GetRadius() << endl;
  cout << "THETA IN: " << theta_in << endl;
  cout << "PHI IN: " << phi_in << endl;

  ev.GetP3(p3);
  ev.GetP4(p4);
  ev_approx.GetP3(p3_approx);
  ev_approx.GetP4(p4_approx);

  Double_t x3_t = p3[0];
  Double_t y3_t = p3[1];
  Double_t x4_t = p4[0];
  Double_t y4_t = p4[1];

  Double_t x3_approx = p3_approx[0];
  Double_t x4_approx = p4_approx[0];
  Double_t y3_approx = p3_approx[1];
  Double_t y4_approx = p4_approx[1];

  cout << "x3_t: " << x3_t << endl;
  cout << "x3_p: " << x3_approx << endl;

  cout << "x4_t: " << x4_t << endl;
  cout << "x4_p: " << x4_approx << endl;

  cout << "y3_t: " << y3_t << endl;
  cout << "y3_p: " << y3_approx << endl;

  cout << "y4_t: " << y4_t << endl;
  cout << "y4_p: " << y4_approx << endl;

  cout << "theta T/Approx: " << theta_out << '\t' << ev_approx.GetTheta_out() << endl;
  cout << "phi T/Approx: " << phi_out << '\t' << ev_approx.GetPhi_out() << endl;
 
  
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

  reco = charge*helix::B*B_event::L/reco;
  reco = reco*0.299792458; //UM

  reco2 = 0.299792458*charge*helix::B*B_event::L/(Sqrt(Sin(theta_out)*Sin(theta_out)-Sin(theta_in)*Sin(theta_in)*Cos(phi_in)*Cos(phi_in))*Sign(1.0,charge*helix::B)-Sin(theta_in)*Sin(phi_in));


  
  cout << "p = " << reco << '\t' << reco2 << endl;
  cout << "p_t = " << p << endl;

}
