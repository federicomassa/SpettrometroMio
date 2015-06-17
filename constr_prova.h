#ifndef CONSTR_PROVA_H
#define CONSTR_PROVA_H

#include <TMatrixD.h>
#include "B_event_approx.h"

Double_t z1 = 50;
Double_t z2 = 60;


TMatrixD Constr(Double_t* vars, Double_t* pars) {

  TMatrixD c(2,1);

  c(0,0) = vars[0]*(z2 - pars[0]) - vars[2]*(z1 - pars[0]);
  c(1,0) = vars[1]*(z2 - pars[0]) - vars[3]*(z1 - pars[0]);
  
  return c;
}

TMatrixD Der_Matrix(Double_t* vars, Double_t* pars) {

  TMatrixD d(4,2);

  d(0,0) = z2 - pars[0];
  d(1,0) = 0;
  d(2,0) = pars[0] - z1;
  d(3,0) = 0;

  d(0,1) = 0;
  d(1,1) = z2 - pars[0];
  d(2,1) = 0;
  d(3,1) = pars[0] - z1;

  // delete warnings
  d(0,0) += 0*vars[0];

  return d;
}

TMatrixD PDer_Matrix(Double_t* vars, Double_t* pars) {

  TMatrixD p(1,2);

  p(0,0) = vars[2] - vars[0];
  p(0,1) = vars[3] - vars[1];

  //delete warnings
  p(0,0) += 0*pars[0];

  
  return p;
}

#endif
