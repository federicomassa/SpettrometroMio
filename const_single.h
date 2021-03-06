#ifndef CONSTR_SINGLE_H
#define CONSTR_SINGLE_H

#include <TMath.h>
#include <TMatrixD.h>
#include "B_event_approx.h"

#ifndef ZARRAY
#define ZARRAY
static Double_t z1 = 50;
static Double_t z2 = 60;
static Double_t z3 = z2 + B_event_approx::L;
static Double_t z4 = z3 + B_event_approx::Delta_z;
#endif

/*
  Vars: x1, y1, x2, y2.... x4, y4
  Pars: z_v, q/p
 */

TMatrixD Constr(Double_t* var, Double_t* par) {

  TMatrixD c(6,1);
  
  c(0,0) = var[0]*(z2 - par[0]) - var[2]*(z1 - par[0]);
  c(1,0) = var[0]*(z3 - par[0]) - var[4]*(z1 - par[0]);
  c(2,0) = var[0]*(z4 - par[0]) - var[6]*(z1 - par[0]);
  c(3,0) = var[1]*(z2 - par[0]) - var[3]*(z1 - par[0]);
  c(4,0) = var[5] - var[1] - 0.5*(var[3]-var[1])*(z2+z3-2*z1)/(z2-z1) - 0.5*(z3-z2)*((var[3]-var[1])/(z2-z1)+B_event_approx::p_kick*par[1]);
  c(5,0) = (var[7]-var[5])/(z4-z3) - (var[3]-var[1])/(z2-z1) - B_event_approx::p_kick*par[1];
  
  return c;

}

TMatrixD Der(Double_t* var, Double_t* par) {

  TMatrixD B(8,6);
  
  for (UInt_t i = 0; i < 8; i++) {
    for (UInt_t j = 0; j < 6; j++) {
      B(i,j) = 0.0*var[i]*par[j]; //delete warnings
    }
  }

  B(0,0) = z2 - par[0];
  B(2,0) = -(z1 - par[0]);
  B(0,1) = z3 - par[0];
  B(4,1) = -(z1 - par[0]);
  B(0,2) = z4 - par[0];
  B(6,2) = -(z1 - par[0]);
  B(1,3) = z2 - par[0];
  B(3,3) = -(z1 - par[0]);
  B(1,4) = -1 + (z3-z1)/(z2-z1);
  B(3,4) = -(z3-z1)/(z2-z1);
  B(5,4) = 1;
  B(1,5) = 1/(z2-z1);
  B(3,5) = -1/(z2-z1);
  B(5,5) = -1/(z4-z3);
  B(7,5) = 1/(z4-z3);
  
  return B;

}

TMatrixD PDer(Double_t* var, Double_t* par) {

  TMatrixD A(2,6);

  for (UInt_t i = 0; i < 2; i++) {
    for (UInt_t j = 0; j < 6; j++) {
      A(i,j) = 0.0*par[i]; //delete warnings
    }
  }

  A(0,0) = var[2]-var[0];
  A(0,1) = var[4]-var[0];
  A(0,2) = var[6]-var[0];
  A(0,3) = var[3]-var[1];
  A(1,4) = -0.5*B_event_approx::p_kick*(z3-z2);
  A(1,5) = -B_event_approx::p_kick;

  return A;

}


#endif
