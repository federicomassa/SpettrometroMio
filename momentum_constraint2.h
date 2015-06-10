#ifndef MOMENTUM_CONSTRAINT_H
#define MOMENTUM_CONSTRAINT_H

  /*
    Variables order
    0: px_plus
    1: py_plus
    2: pz_plus
    3: px_min
    4: py_min
    5: pz_min
   */

#include <TMatrixD.h>
#include <TMath.h>

using namespace TMath;


Double_t Energy(Double_t mass, Double_t px, Double_t py, Double_t pz) {
  return Sqrt(px*px + py*py + pz*pz + mass*mass);
}

TMatrixD Momentum_Constr_Vector(Double_t* p) {
 
  Double_t pi_mass = 0.13957018; // Charged Pi Mass
  
  TMatrixD p_Vector(2,1);
  
  p_Vector(0,0) = p[0] + p[3];
  p_Vector(1,0) = p[1] + p[4];

  return p_Vector;
}

TMatrixD Momentum_Constr_Derivative(Double_t* p) {
  
  Double_t pi_mass = 0.13957018;

  TMatrixD D(6,2);
  
  for (UInt_t i = 0; i < 6; i++) {
    for (UInt_t j = 0; j < 2; j++) {
      D(i,j) = 0.0;
    }
  }

  D(0,0) = 1;
  D(3,0) = 1;
  D(1,1) = 1;
  D(4,1) = 1;

  return D;
}



#endif
