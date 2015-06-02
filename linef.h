#ifndef LINE_H
#define LINE_H

#include "linef.cpp"
#include <TMath.h>

line::line() {
}

line::line(Float_t* p1, Float_t* p2) {

 Float_t norm2;

  norm2 = TMath::Power(p2[0]-p1[0],2) + TMath::Power(p2[1]-p1[1],2) + TMath::Power(p2[2]-p1[2],2);

  init_point[0] = p1[0];
  init_point[1] = p1[1];
  init_point[2] = p1[2];

  

  theta = TMath::ACos((p2[2]-p1[2])/TMath::Sqrt(norm2));

  phi = TMath::ATan((p2[1]-p1[1])/(p2[0]-p1[0]));
  // ATan gives -pi/2 < phi < pi/2
  if (phi < 0) phi += 2*TMath::Pi();
  if (Float_t(TMath::Sign(1.0,p2[0]-p1[0]))*Float_t(TMath::Sign(1.0,TMath::Sin(theta))) < 0) phi += TMath::Pi();

  coeff[0] = TMath::Sin(theta)*TMath::Cos(phi);
  coeff[1] = TMath::Sin(theta)*TMath::Sin(phi);
  coeff[2] = TMath::Cos(theta);
}

line::line(Float_t* p, Float_t theta0, Float_t phi0) {
  init_point[0] = p[0];
  init_point[1] = p[1];
  init_point[2] = p[2];

  theta = theta0;
  phi = phi0;

  coeff[0] = TMath::Sin(theta0)*TMath::Cos(phi0);
  coeff[1] = TMath::Sin(theta0)*TMath::Sin(phi0);
  coeff[2] = TMath::Cos(theta0);
}

void line::SetLine(Float_t* p1, Float_t* p2) {
  
  Float_t norm2;

  norm2 = TMath::Power(p2[0]-p1[0],2) + TMath::Power(p2[1]-p1[1],2) + TMath::Power(p2[2]-p1[2],2);

  init_point[0] = p1[0];
  init_point[1] = p1[1];
  init_point[2] = p1[2];

  

  theta = TMath::ACos((p2[2]-p1[2])/TMath::Sqrt(norm2));

  phi = TMath::ATan((p2[1]-p1[1])/(p2[0]-p1[0]));
  // ATan gives -pi/2 < phi < pi/2
  if (phi < 0) phi += 2*TMath::Pi();
  if (Float_t(TMath::Sign(1.0,p2[0]-p1[0]))*Float_t(TMath::Sign(1.0,TMath::Sin(theta))) < 0) phi += TMath::Pi();

  coeff[0] = TMath::Sin(theta)*TMath::Cos(phi);
  coeff[1] = TMath::Sin(theta)*TMath::Sin(phi);
  coeff[2] = TMath::Cos(theta);
}

void line::SetLine(Float_t* p, Float_t theta0, Float_t phi0) {

  init_point[0] = p[0];
  init_point[1] = p[1];
  init_point[2] = p[2];

  theta = theta0;
  phi = phi0;

  coeff[0] = TMath::Sin(theta0)*TMath::Cos(phi0);
  coeff[1] = TMath::Sin(theta0)*TMath::Sin(phi0);
  coeff[2] = TMath::Cos(theta0);
}

//Gives t
Float_t line::GetParameterAtZ(Float_t z) {
  return (z-init_point[2])/coeff[2];
}

void line::GetPoint(Float_t t, Float_t* p)  {
  p[0] = init_point[0] + coeff[0]*t;
  p[1] = init_point[1] + coeff[1]*t;
  p[2] = init_point[2] + coeff[2]*t;
}

void line::GetInitCoords(Float_t* p0) {
  p0[0] = init_point[0];
  p0[1] = init_point[1];
  p0[2] = init_point[2];
}

void line::GetCoefficients(Float_t* coeff0) {
  coeff0[0] = coeff[0];
  coeff0[1] = coeff[1];
  coeff0[2] = coeff[2];
}

Float_t line::GetTheta() {
  return theta;
}

Float_t line::GetPhi() {
  return phi;
}

#endif
