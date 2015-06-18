#ifndef CONSTR_DOUBLE_H
#define CONSTR_DOUBLE_H

#include <TMath.h>
#include <TMatrixD.h>
#include "B_event_approx.h"


#ifndef ZARRAY
#define ZARRAY
Double_t z1 = 50;
Double_t z2 = 60;
Double_t z3 = z2 + B_event_approx::L;
Double_t z4 = z3 + B_event_approx::Delta_z;
#endif

#ifndef MASSES
#define MASSES
Double_t K_Mass = 0.497614;
Double_t Pi_Mass = 0.13957018;
#endif


/*
  Vars: x1, y1, x2, y2.... x4, y4 for + and -
  Pars: z_v, q1/p, q2/p
 */

Double_t Delta_x_plus(Double_t* var) {
  return (var[2]-var[0]);
}

Double_t Delta_y_plus(Double_t* var) {
  return (var[3]-var[1]);
}

Double_t Delta_x_min(Double_t* var) {
  return (var[10]-var[8]);
}

Double_t Delta_y_min(Double_t* var) {
  return (var[11]-var[9]);
}

Double_t Norm_plus(Double_t* var) {
  return TMath::Sqrt(TMath::Power(Delta_x_plus(var),2) + TMath::Power(Delta_y_plus(var),2) + TMath::Power(z2-z1,2));
}

Double_t Norm_min(Double_t* var) {
  return TMath::Sqrt(TMath::Power(Delta_x_min(var),2) + TMath::Power(Delta_y_min(var),2) + TMath::Power(z2-z1,2));
}

Double_t Energy_plus(Double_t* par) {
  return TMath::Sqrt(1/par[1]/par[1] + Pi_Mass*Pi_Mass);
}

Double_t Energy_min(Double_t* par) {
  return TMath::Sqrt(1/par[2]/par[2] + Pi_Mass*Pi_Mass);
}

Double_t Invariant_Mass(Double_t* var, Double_t* par) {

  Double_t W;
  
  W = Pi_Mass*Pi_Mass - K_Mass*K_Mass/2 + TMath::Sqrt(1/(par[1]*par[1]) + Pi_Mass*Pi_Mass)*TMath::Sqrt(1/(par[2]*par[2]) + Pi_Mass*Pi_Mass) + 1/(par[1]*par[1])*(Norm_plus(var)*Norm_plus(var)- (z2-z1)*(z2-z1))/(Norm_plus(var)*Norm_plus(var)) + 1/(par[1]*par[2])*TMath::Power(z2-z1,2)/Norm_plus(var)/Norm_min(var);

  return W;
}

TMatrixD Constr(Double_t* var, Double_t* par) {

  TMatrixD c(15,1);
  
  c(0,0) = var[0]*(z2 - par[0]) - var[2]*(z1 - par[0]);
  c(1,0) = var[0]*(z3 - par[0]) - var[4]*(z1 - par[0]);
  c(2,0) = var[0]*(z4 - par[0]) - var[6]*(z1 - par[0]);
  c(3,0) = var[1]*(z2 - par[0]) - var[3]*(z1 - par[0]);
  c(4,0) = var[5] - var[1] - 0.5*(var[3]-var[1])*(z2+z3-2*z1)/(z2-z1) - 0.5*(z3-z2)*((var[3]-var[1])/(z2-z1)+B_event_approx::p_kick*par[1]);
  c(5,0) = (var[7]-var[5])/(z4-z3) - (var[3]-var[1])/(z2-z1) - B_event_approx::p_kick*par[1];

  c(6,0) = var[8]*(z2 - par[0]) - var[10]*(z1 - par[0]);
  c(7,0) = var[8]*(z3 - par[0]) - var[12]*(z1 - par[0]);
  c(8,0) = var[8]*(z4 - par[0]) - var[14]*(z1 - par[0]);
  c(9,0) = var[9]*(z2 - par[0]) - var[11]*(z1 - par[0]);
  c(10,0) = var[13] - var[9] - 0.5*(var[11]-var[9])*(z2+z3-2*z1)/(z2-z1) - 0.5*(z3-z2)*((var[11]-var[9])/(z2-z1)+B_event_approx::p_kick*par[2]);
  c(11,0) = (var[15]-var[13])/(z4-z3) - (var[11]-var[9])/(z2-z1) - B_event_approx::p_kick*par[2];

  c(12,0) = par[2]*Delta_x_plus(var)*Norm_min(var) - par[1]*Delta_x_min(var)*Norm_plus(var);
  c(13,0) = par[2]*Delta_y_plus(var)*Norm_min(var) - par[1]*Delta_y_min(var)*Norm_plus(var);

  c(14,0) = Invariant_Mass(var, par);
 
  return c;

}

TMatrixD Der(Double_t* var, Double_t* par) {

  TMatrixD B(16,15);
  
  for (UInt_t i = 0; i < 16; i++) {
    for (UInt_t j = 0; j < 15; j++) {
      B(i,j) = 0.0*var[i]; //delete warnings
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

  B(8,6) = z2 - par[0];
  B(10,6) = -(z1 - par[0]);
  B(8,7) = z3 - par[0];
  B(12,7) = -(z1 - par[0]);
  B(8,8) = z4 - par[0];
  B(14,8) = -(z1 - par[0]);
  B(9,9) = z2 - par[0];
  B(11,9) = -(z1 - par[0]);
  B(9,10) = -1 + (z3-z1)/(z2-z1);
  B(11,10) = -(z3-z1)/(z2-z1);
  B(13,10) = 1;
  B(9,11) = 1/(z2-z1);
  B(11,11) = -1/(z2-z1);
  B(13,11) = -1/(z4-z3);
  B(15,11) = 1/(z4-z3);

  B(0,12) = -par[2]*Norm_min(var) + par[1]*Delta_x_min(var)*Delta_x_plus(var)/Norm_plus(var);
  B(1,12) = par[1]*Delta_x_min(var)*Delta_y_plus(var)/Norm_plus(var);
  B(2,12) = -B(0,12);
  B(3,12) = -B(1,12);
  B(8,12) = par[1]*Norm_plus(var) - par[2]*Delta_x_min(var)*Delta_x_plus(var)/Norm_min(var);
  B(9,12) = -par[2]*Delta_x_plus(var)*Delta_y_min(var)/Norm_min(var);
  B(10,12) = -B(8,12);
  B(11,12) = -B(9,12);


  B(0,13) = par[1]*Delta_y_min(var)*Delta_x_plus(var)/Norm_plus(var);
  B(1,13) = -par[2]*Norm_min(var) + par[1]*Delta_y_min(var)*Delta_y_plus(var)/Norm_plus(var);
  B(2,13) = -B(0,13);
  B(3,13) = -B(1,13);
  B(8,13) = -par[2]*Delta_y_plus(var)*Delta_x_min(var)/Norm_min(var);
  B(9,13) = par[1]*Norm_plus(var) - par[2]*Delta_y_min(var)*Delta_y_plus(var)/Norm_min(var);
  B(10,13) = -B(8,13);
  B(11,13) = -B(9,13);

  B(0,14) = -1/(par[1]*par[1])*2*TMath::Power(z2-z1,2)*Delta_x_plus(var)/TMath::Power(Norm_plus(var),4) + 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_x_plus(var)/Norm_min(var)/TMath::Power(Norm_plus(var),3);
  B(1,14) = -1/(par[1]*par[1])*2*TMath::Power(z2-z1,2)*Delta_y_plus(var)/TMath::Power(Norm_plus(var),4) + 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_y_plus(var)/Norm_min(var)/TMath::Power(Norm_plus(var),3);
  B(2,14) = -B(0,14);
  B(3,14) = -B(1,14);
  B(8,14) = 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_x_min(var)/Norm_plus(var)/TMath::Power(Norm_min(var),3);
  B(9,14) = 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_y_min(var)/Norm_plus(var)/TMath::Power(Norm_min(var),3);
  B(10,14) = -B(8,14);
  B(11,14) = -B(9,14);

  
  return B;

}

TMatrixD PDer(Double_t* var, Double_t* par) {

  TMatrixD A(3,15);

  for (UInt_t i = 0; i < 3; i++) {
    for (UInt_t j = 0; j < 15; j++) {
      A(i,j) = 0.0*par[i]; //delete warnings
    }
  }

  A(0,0) = var[2]-var[0];
  A(0,1) = var[4]-var[0];
  A(0,2) = var[6]-var[0];
  A(0,3) = var[3]-var[1];

  A(0,6) = var[10]-var[8];
  A(0,7) = var[12]-var[8];
  A(0,8) = var[14]-var[8];
  A(0,9) = var[11]-var[9];

  A(1,4) = -0.5*B_event_approx::p_kick*(z3-z2);
  A(1,5) = -B_event_approx::p_kick;
  A(1,12) = -Delta_x_min(var)*Norm_plus(var);
  A(1,13) = -Delta_y_min(var)*Norm_plus(var);;
  A(1,14) = -Energy_min(par)/Energy_plus(par)*TMath::Power(1/par[1],3) - 2*TMath::Power(1/par[1],3)*(TMath::Power(Norm_plus(var),2)-TMath::Power(z2-z1,2))-TMath::Power(1/par[1],2)/par[2]*TMath::Power(z2-z1,2)/Norm_plus(var)/Norm_min(var);

  A(2,10) = -0.5*B_event_approx::p_kick*(z3-z2);
  A(2,11) = -B_event_approx::p_kick;
  A(2,12) = Delta_x_plus(var)*Norm_min(var);
  A(2,13) = Delta_y_plus(var)*Norm_min(var);;
  A(2,14) = -Energy_plus(par)/Energy_min(par)*TMath::Power(1/par[2],3) - TMath::Power(1/par[2],2)/par[1]*TMath::Power(z2-z1,2)/Norm_plus(var)/Norm_min(var);

  return A;

}


#endif
