#ifndef GETTHETAPHI_H
#define GETTHETAPHI_H

#include <TMath.h>

using namespace TMath;

void GetThetaPhi(Double_t* p1, Double_t* p2, Double_t &theta, Double_t &phi) {

  Double_t norm, norm2;

  norm2 = Power(p2[0]-p1[0],2) + Power(p2[1]-p1[1],2) + Power(p2[2]-p1[2],2);
  norm = Sqrt(norm2);

  theta = ACos((p2[2]-p1[2])/norm);
  phi = ATan((p2[1]-p1[1])/(p2[0]-p1[0]));

  if (Sign(1.0, p2[0]-p1[0])*Sign(1.0,Sin(theta)) < 0) phi += Pi();
  if (phi < 0) phi += 2*Pi();
  
}

#endif
