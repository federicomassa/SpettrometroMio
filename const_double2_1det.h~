#ifndef CONSTR_DOUBLE2_1DET_H
#define CONSTR_DOUBLE2_1DET_H

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
  Vars: x1, y1, x2, y2, x3, y3, for + and -
  Pars: z_v, q1/p, q2/p
 */

Double_t Delta_x_plus(Double_t* var) {
  return (var[2]-var[0]);
}

Double_t Delta_y_plus(Double_t* var) {
  return (var[3]-var[1]);
}

Double_t Delta_x_min(Double_t* var) {
  return (var[8]-var[6]);
}

Double_t Delta_y_min(Double_t* var) {
  return (var[9]-var[7]);
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

  TMatrixD c(11,1);
  
  c(0,0) = var[0]*(z2 - par[0]) - var[2]*(z1 - par[0]);
  c(1,0) = var[0]*(z3 - par[0]) - var[4]*(z1 - par[0]);
  c(2,0) = var[1]*(z2 - par[0]) - var[3]*(z1 - par[0]);
  c(3,0) = var[5] - var[1] - 0.5*(var[3]-var[1])*(z2+z3-2*z1)/(z2-z1) - 0.5*(z3-z2)*((var[3]-var[1])/(z2-z1)+B_event_approx::p_kick*par[1]);


  c(4,0) = var[6]*(z2 - par[0]) - var[8]*(z1 - par[0]);
  c(5,0) = var[6]*(z3 - par[0]) - var[10]*(z1 - par[0]);
  c(6,0) = var[7]*(z2 - par[0]) - var[9]*(z1 - par[0]);
  c(7,0) = var[11] - var[7] - 0.5*(var[9]-var[7])*(z2+z3-2*z1)/(z2-z1) - 0.5*(z3-z2)*((var[9]-var[7])/(z2-z1)+B_event_approx::p_kick*par[2]);


  c(8,0) = par[2]*Delta_x_plus(var)*Norm_min(var) - par[1]*Delta_x_min(var)*Norm_plus(var);
  c(9,0) = par[2]*Delta_y_plus(var)*Norm_min(var) - par[1]*Delta_y_min(var)*Norm_plus(var);

  c(10,0) = Invariant_Mass(var, par);
 
  return c;

}

TMatrixD Der(Double_t* var, Double_t* par) {

  TMatrixD B(12,11);
  
  for (UInt_t i = 0; i < 12; i++) {
    for (UInt_t j = 0; j < 11; j++) {
      B(i,j) = 0.0*var[i]; //delete warnings
    }
  }

  B(0,0) = z2 - par[0];
  B(2,0) = -(z1 - par[0]);
  B(0,1) = z3 - par[0];
  B(4,1) = -(z1 - par[0]);
  B(1,2) = z2 - par[0];
  B(3,2) = -(z1 - par[0]);
  B(1,3) = -1 + (z3-z1)/(z2-z1);
  B(3,3) = -(z3-z1)/(z2-z1);
  B(5,3) = 1;


  B(6,4) = z2 - par[0];
  B(8,4) = -(z1 - par[0]);
  B(6,5) = z3 - par[0];
  B(10,5) = -(z1 - par[0]);
  B(7,6) = z2 - par[0];
  B(9,6) = -(z1 - par[0]);
  B(7,7) = -1 + (z3-z1)/(z2-z1);
  B(9,7) = -(z3-z1)/(z2-z1);
  B(11,7) = 1;

  B(0,8) = -par[2]*Norm_min(var) + par[1]*Delta_x_min(var)*Delta_x_plus(var)/Norm_plus(var);
  B(1,8) = par[1]*Delta_x_min(var)*Delta_y_plus(var)/Norm_plus(var);
  B(2,8) = -B(0,8);
  B(3,8) = -B(1,8);
  B(6,8) = par[1]*Norm_plus(var) - par[2]*Delta_x_min(var)*Delta_x_plus(var)/Norm_min(var);
  B(7,8) = -par[2]*Delta_x_plus(var)*Delta_y_min(var)/Norm_min(var);
  B(8,8) = -B(6,8);
  B(9,8) = -B(7,8);


  B(0,9) = par[1]*Delta_y_min(var)*Delta_x_plus(var)/Norm_plus(var);
  B(1,9) = -par[2]*Norm_min(var) + par[1]*Delta_y_min(var)*Delta_y_plus(var)/Norm_plus(var);
  B(2,9) = -B(0,9);
  B(3,9) = -B(1,9);
  B(6,9) = -par[2]*Delta_y_plus(var)*Delta_x_min(var)/Norm_min(var);
  B(7,9) = par[1]*Norm_plus(var) - par[2]*Delta_y_min(var)*Delta_y_plus(var)/Norm_min(var);
  B(8,9) = -B(6,9);
  B(9,9) = -B(7,9);

  B(0,10) = -1/(par[1]*par[1])*2*TMath::Power(z2-z1,2)*Delta_x_plus(var)/TMath::Power(Norm_plus(var),4) + 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_x_plus(var)/Norm_min(var)/TMath::Power(Norm_plus(var),3);
  B(1,10) = -1/(par[1]*par[1])*2*TMath::Power(z2-z1,2)*Delta_y_plus(var)/TMath::Power(Norm_plus(var),4) + 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_y_plus(var)/Norm_min(var)/TMath::Power(Norm_plus(var),3);
  B(2,10) = -B(0,10);
  B(3,10) = -B(1,10);
  B(6,10) = 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_x_min(var)/Norm_plus(var)/TMath::Power(Norm_min(var),3);
  B(7,10) = 1/(par[1]*par[2])*TMath::Power(z2-z1,2)*Delta_y_min(var)/Norm_plus(var)/TMath::Power(Norm_min(var),3);
  B(8,10) = -B(6,10);
  B(9,10) = -B(7,10);

  
  return B;

}

TMatrixD PDer(Double_t* var, Double_t* par) {

  TMatrixD A(3,11);

  for (UInt_t i = 0; i < 3; i++) {
    for (UInt_t j = 0; j < 11; j++) {
      A(i,j) = 0.0*par[i]; //delete warnings
    }
  }

  A(0,0) = var[2]-var[0];
  A(0,1) = var[4]-var[0];
  A(0,2) = var[3]-var[1];

  A(0,4) = var[8]-var[6];
  A(0,5) = var[10]-var[6];
  A(0,6) = var[9]-var[7];

  A(1,3) = -0.5*B_event_approx::p_kick*(z3-z2);
 
  A(1,8) = -Delta_x_min(var)*Norm_plus(var);
  A(1,9) = -Delta_y_min(var)*Norm_plus(var);;
  A(1,10) = -Energy_min(par)/Energy_plus(par)*TMath::Power(1/par[1],3) - 2*TMath::Power(1/par[1],3)*(TMath::Power(Norm_plus(var),2)-TMath::Power(z2-z1,2))-TMath::Power(1/par[1],2)/par[2]*TMath::Power(z2-z1,2)/Norm_plus(var)/Norm_min(var);

  A(2,7) = -0.5*B_event_approx::p_kick*(z3-z2);
 
  A(2,8) = Delta_x_plus(var)*Norm_min(var);
  A(2,9) = Delta_y_plus(var)*Norm_min(var);;
  A(2,10) = -Energy_plus(par)/Energy_min(par)*TMath::Power(1/par[2],3) - TMath::Power(1/par[2],2)/par[1]*TMath::Power(z2-z1,2)/Norm_plus(var)/Norm_min(var);

  return A;

}


#endif
