// as InitializeZPP but without 4th detector

#ifndef INITIALIZEZPP_1DET_H
#define INITIALIZEZPP_1DET_H

//16 vars.. x1,y1 ... x4,y4  *  2
//3 pars.. z_decay, q+/p+, q-/p-

#include "B_event_approx.h"

#ifndef ZARRAY
#define ZARRAY
Double_t z1 = 50;
Double_t z2 = 60;
Double_t z3 = z2 + B_event_approx::L;
Double_t z4 = z3 + B_event_approx::Delta_z;
#endif


void InitializeZPP_1det(Double_t* init_m, Double_t* init_p) {

  // raw esteem for z: CDA_plus
  init_p[0] = z1 - (z2 - z1)*(init_m[0]*(init_m[2]-init_m[0]) + init_m[1]*(init_m[3]-init_m[1]))/((init_m[2] - init_m[0])*(init_m[2] - init_m[0]) + (init_m[3] - init_m[1])*(init_m[3] - init_m[1]));
  init_p[1] = ((init_m[5]-init_m[1]) - (init_m[3]-init_m[1])*(z3-z1)/(z2-z1))*2/(z3-z2)/B_event_approx::p_kick;
  init_p[2] =  ((init_m[11]-init_m[7]) - (init_m[9]-init_m[7])*(z3-z1)/(z2-z1))*2/(z3-z2)/B_event_approx::p_kick;

  /* if (init_p[1]*init_p[2] > 0) { */
  /*   if (TMath::Abs(init_p[1]) < TMath::Abs(init_p[2])) init_p[1] = -init_p[1]; */
  /*   else init_p[2] = -init_p[2]; */
  /* } */
  
}

#endif
