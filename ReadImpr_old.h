#ifndef READIMPR_H
#define READIMPR_H

#include "GetElement.h"
#include <fstream>
#include <string>

void ReadImpr(ifstream &in, Int_t &ev_no, Double_t &K_z, Double_t &K_p, Double_t &pi_plus_modp, Double_t &pi_plus_theta, Double_t &pi_plus_phi, Double_t &pi_min_modp, Double_t &pi_min_theta, Double_t &pi_min_phi, Double_t* P1_plus, Double_t* P1_min, Double_t* P2_plus, Double_t* P2_min) {

  string str;
  getline(in,str);
  
  ev_no = Int_t(GetElement(str, 0));
  K_z = GetElement(str, 1);
  K_p = GetElement(str, 2);
  pi_plus_modp = GetElement(str, 3);
  pi_plus_theta = GetElement(str, 4);
  pi_plus_phi = GetElement(str, 5);
  pi_min_modp = GetElement(str, 6);
  pi_min_theta = GetElement(str, 7);
  pi_min_phi = GetElement(str,8);
  P1_plus[0] = GetElement(str, 9);
  P1_plus[1] = GetElement(str, 10);
  P1_min[0] = GetElement(str, 11);
  P1_min[1] = GetElement(str, 12);
  P2_plus[0] = GetElement(str, 13);
  P2_plus[1] = GetElement(str, 14);
  P2_min[0] = GetElement(str, 15);
  P2_min[1] = GetElement(str, 16);  
}

#endif
