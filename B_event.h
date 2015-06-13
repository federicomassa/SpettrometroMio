#ifndef B_EVENT_H
#define B_EVENT_H

#include "B_event.cpp"
#include "line.h"
#include "helix.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TCanvas.h>

B_event::B_event() {}

B_event::B_event(Double_t* p10, Double_t* p20, Double_t charge0, Double_t gamma0) {

  Double_t vector_out[3]; //outgoing particle versor
  Double_t beta;

  p1[0] = p10[0];
  p1[1] = p10[1];
  p1[2] = p10[2];

  p2[0] = p20[0];
  p2[1] = p20[1];
  p2[2] = p20[2];

  charge = charge0;
  gamma = gamma0;

  beta = TMath::Sqrt(gamma0*gamma0-1)/gamma; 

  line_in.SetLine(p10, p20);

  theta_in = line_in.GetTheta();
  phi_in = line_in.GetPhi();

  /* //APPROXIMATION: ultrarelativistic particle */
  /* init_speed[0] = TMath::Sin(theta_in)*TMath::Cos(phi_in);  */
  /* init_speed[1] = TMath::Sin(theta_in)*TMath::Sin(phi_in); */
  /* init_speed[2] = TMath::Cos(theta_in); */

  init_speed[0] = beta*TMath::Sin(theta_in)*TMath::Cos(phi_in);
  init_speed[1] = beta*TMath::Sin(theta_in)*TMath::Sin(phi_in);
  init_speed[2] = beta*TMath::Cos(theta_in);
  

  h.SetHelix(p20, init_speed, charge0, gamma0);

  omega = h.GetOmega();
  radius = h.GetRadius();
  
  cout << "OMEGA: " << omega << endl;

  /* time_interval = L/TMath::C()*1/TMath::Sqrt(1-TMath::Power(TMath::Sin(theta_in)*TMath::Cos(phi_in),2)); //APPROXIMATION at small theta */

  omegaT = (TMath::ASin((L/radius)*TMath::Sign(1.0,charge0*helix::B) + TMath::Sin(theta_in)*TMath::Sin(phi_in)/TMath::Sqrt(1-TMath::Power(TMath::Sin(theta_in)*TMath::Cos(phi_in),2))) - TMath::ATan(TMath::Tan(theta_in)*TMath::Sin(phi_in)));

  time_interval = omegaT/omega;

  h.GetPoint(time_interval, p3);
  
  vector_out[0] = init_speed[0]/beta ;
  vector_out[1] = (init_speed[1]*TMath::Cos(omegaT) + init_speed[2]*TMath::Sin(omegaT))/beta;
  vector_out[2] = (init_speed[2]*TMath::Cos(omegaT) - init_speed[1]*TMath::Sin(omegaT))/beta;

  //THETA-PHI OUT//
  theta_out = TMath::ACos(vector_out[2]);
  phi_out = TMath::ATan(vector_out[1]/vector_out[0]);

  if (phi_out < 0) phi_out += 2*TMath::Pi();
  if (TMath::Sign(1.0,vector_out[0])*TMath::Sign(1.0,TMath::Sin(theta_out)) < 0) phi_out += TMath::Pi();
  // // // //

  line_out.SetLine(p3, theta_out, phi_out);  
  
  line_out.GetPoint(line_out.GetParameterAtZ(p2[2]+L+Delta_z), p4);

}

//As constructor
void B_event::SetB_event(Double_t* p10, Double_t* p20, Double_t charge0, Double_t gamma0) {

  Double_t vector_out[3]; //outgoing particle versor
  Double_t beta;

  p1[0] = p10[0];
  p1[1] = p10[1];
  p1[2] = p10[2];

  p2[0] = p20[0];
  p2[1] = p20[1];
  p2[2] = p20[2];

  charge = charge0;
  gamma = gamma0;

  beta = TMath::Sqrt(gamma0*gamma0-1)/gamma; 

  line_in.SetLine(p10, p20);

  theta_in = line_in.GetTheta();
  phi_in = line_in.GetPhi();

  /* //APPROXIMATION: ultrarelativistic particle */
  /* init_speed[0] = TMath::Sin(theta_in)*TMath::Cos(phi_in);  */
  /* init_speed[1] = TMath::Sin(theta_in)*TMath::Sin(phi_in); */
  /* init_speed[2] = TMath::Cos(theta_in); */

  init_speed[0] = beta*TMath::Sin(theta_in)*TMath::Cos(phi_in);
  init_speed[1] = beta*TMath::Sin(theta_in)*TMath::Sin(phi_in);
  init_speed[2] = beta*TMath::Cos(theta_in);
  

  h.SetHelix(p20, init_speed, charge0, gamma0);

  omega = h.GetOmega();
  radius = h.GetRadius();

  /* time_interval = L/TMath::C()*1/TMath::Sqrt(1-TMath::Power(TMath::Sin(theta_in)*TMath::Cos(phi_in),2)); //APPROXIMATION at small theta */

  omegaT = (TMath::ASin((L/radius)*TMath::Sign(1.0,charge0*helix::B) + TMath::Sin(theta_in)*TMath::Sin(phi_in)/TMath::Sqrt(1-TMath::Power(TMath::Sin(theta_in)*TMath::Cos(phi_in),2))) - TMath::ATan(TMath::Tan(theta_in)*TMath::Sin(phi_in)));

  time_interval = omegaT/omega;

  h.GetPoint(time_interval, p3);
  
  vector_out[0] = init_speed[0]/beta ;
  vector_out[1] = (init_speed[1]*TMath::Cos(omegaT) + init_speed[2]*TMath::Sin(omegaT))/beta;
  vector_out[2] = (init_speed[2]*TMath::Cos(omegaT) - init_speed[1]*TMath::Sin(omegaT))/beta;

  //THETA-PHI OUT//
  theta_out = TMath::ACos(vector_out[2]);
  phi_out = TMath::ATan(vector_out[1]/vector_out[0]);

  if (phi_out < 0) phi_out += 2*TMath::Pi();
  if (TMath::Sign(1.0,vector_out[0])*TMath::Sign(1.0,TMath::Sin(theta_out)) < 0) phi_out += TMath::Pi();
  // // // //

  line_out.SetLine(p3, theta_out, phi_out);  
  
  line_out.GetPoint(line_out.GetParameterAtZ(p2[2]+L+Delta_z), p4);

}

TGraph2D* B_event::GetEventGraph() {

  Double_t par_in_max, par_out_max;
  Double_t point_buff[3];
  Double_t final_point[3];
  Double_t line_out_length = 30;
  const Int_t npoints = 10000;
  
  line_out.GetPoint(line_out_length, final_point);

  TGraph2D* event_display = new TGraph2D(3*(npoints+1));
  event_display->SetTitle("Event Display;x;y;z");

  par_in_max = line_in.GetParameterAtZ(p2[2]);
  par_out_max = line_out.GetParameterAtZ(final_point[2]);
  
  

  // DRAW line_in
  for (int i = 0; i < npoints + 1; i++) {
    line_in.GetPoint(par_in_max/Double_t(npoints)*Double_t(i), point_buff);
    event_display->SetPoint(i+1, point_buff[0], point_buff[1], point_buff[2]);
  }
  
  // DRAW helix
  for (int i = 0; i < npoints + 1; i++) {
    h.GetPoint(time_interval/Double_t(npoints)*Double_t(i), point_buff);
    event_display->SetPoint(npoints + 1 + i+1, point_buff[0], point_buff[1], point_buff[2]);
  }

  /* cout << "LAST HELIX Z: " << point_buff[2] << endl; */

  // DRAW line_out
  for (int i = 0; i < npoints + 1; i++) {
    line_out.GetPoint(par_out_max/Double_t(npoints)*Double_t(i), point_buff);
    event_display->SetPoint(2*(npoints+1) + i+1, point_buff[0], point_buff[1], point_buff[2]);
  }

  /* cout << "FINAL Z: " << final_point[2] << endl; */

  return event_display;

}

Double_t B_event::GetTheta_in() {
  return theta_in;
}

Double_t B_event::GetTheta_out() {
  return theta_out;
}

Double_t B_event::GetPhi_in() {
  return phi_in;
}

Double_t B_event::GetPhi_out() {
  return phi_out;
}

Double_t B_event::GetTimeInterval() {
  return time_interval;
}

Double_t B_event::GetOmega() {
  return omega;
}

Double_t B_event::GetOmegaT() {
  return omegaT;
}

Double_t B_event::GetRadius() {
  return radius;
}

void B_event::GetP3(Double_t* p30) {
  p30[0] = p3[0];
  p30[1] = p3[1];
  p30[2] = p3[2];
}

void B_event::GetP4(Double_t* p40) {
  p40[0] = p4[0];
  p40[1] = p4[1];
  p40[2] = p4[2];
}

Double_t B_event::L = 0.35;
Double_t B_event::Delta_z = 10;



#endif
