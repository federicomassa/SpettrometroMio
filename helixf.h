#ifndef HELIX_H
#define HELIX_H

#include "helixf.cpp"
#include <TMath.h>
#include <iostream>

helix::helix() {};

//charge in e units, coordinates in meters, speed in units of c
helix::helix(Float_t* init_p, Float_t* init_sp, Float_t charge0, Float_t gamma0) {
  
  init_point[0] = init_p[0];
  init_point[1] = init_p[1];
  init_point[2] = init_p[2];
  
  init_speed[0] = init_sp[0];
  init_speed[1] = init_sp[1];
  init_speed[2] = init_sp[2];

  charge = charge0;

  gamma = gamma0;

  // [omega] = s^-1
  omega = TMath::C()*TMath::C()/TMath::Power(10,9)*charge0*B/(mass*gamma0);
}


void helix::SetHelix(Float_t* init_p, Float_t* init_sp, Float_t charge0, Float_t gamma0) {

  init_point[0] = init_p[0];
  init_point[1] = init_p[1];
  init_point[2] = init_p[2];

  init_speed[0] = init_sp[0];
  init_speed[1] = init_sp[1];
  init_speed[2] = init_sp[2];

  charge = charge0;

  gamma = gamma0;

  // [omega] = s^-1
  omega = TMath::C()*TMath::C()/TMath::Power(10,9)*charge0*B/(mass*gamma0);

}

void helix::GetInitCoord(Float_t *init_p) {
  init_p[0] = init_point[0];
  init_p[1] = init_point[1];
  init_p[2] = init_point[2];
}

void helix::GetInitSpeed(Float_t *init_sp) {
  init_sp[0] = init_speed[0];
  init_sp[1] = init_speed[1];
  init_sp[2] = init_speed[2];
}

Float_t helix::GetCharge() {
  return charge;
}

Float_t helix::GetGamma() {
  return gamma;
}

Float_t helix::GetOmega() {
  return omega;
}

Float_t helix::GetRadius() {
  return TMath::Abs(1/omega*TMath::C()*TMath::Sqrt(init_speed[1]*init_speed[1] + init_speed[2]*init_speed[2]));
}


Float_t helix::GetParameterAtX(Float_t x) {

  if (init_speed[0] == 0) std::cout << "ERROR in GetParameter, x_speed = 0" << std::endl;

  return (x-init_point[0])/init_speed[0];
}



void helix::GetPoint(Float_t par, Float_t* p) {
  

  p[0] = init_point[0] + TMath::C()*init_speed[0]*par;
  
  p[1] = init_point[1] + TMath::C()/omega*(init_speed[2] +
				  init_speed[1]*TMath::Sin(omega*par) - 
				  init_speed[2]*TMath::Cos(omega*par));

  p[2] = init_point[2] + TMath::C()/omega*(-init_speed[1] +
				  init_speed[2]*TMath::Sin(omega*par) + 
				  init_speed[1]*TMath::Cos(omega*par));

}

////////// EXPERIMENTAL SETUP //////////////////
Float_t helix::mass = 0.13957; // [GeV/c^2]
Float_t helix::B = 1; // [T]

#endif
