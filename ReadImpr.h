#ifndef READIMPR_H
#define READIMPR_H

#include "GetElement.h"
#include <TNtupleD.h>
#include <string>

void ReadImpr(TNtupleD* nt, Int_t entry, Double_t &K_z, Double_t &K_p, Double_t &pi_plus_modp, Double_t &pi_plus_theta, Double_t &pi_plus_phi, Double_t &pi_min_modp, Double_t &pi_min_theta, Double_t &pi_min_phi, Double_t* P1_plus, Double_t* P1_min, Double_t* P2_plus, Double_t* P2_min) {
  
  Double_t* args;
  nt->GetEntry(entry);
  args = (Double_t*) nt->GetArgs();
  K_z = args[0];
  K_p = args[1];
  pi_plus_modp = args[2];
  pi_plus_theta = args[3];
  pi_plus_phi = args[4];
  pi_min_modp = args[5];
  pi_min_theta = args[6];
  pi_min_phi = args[7];
  P1_plus[0] = args[8];
  P1_plus[1] = args[9];
  P1_min[0] = args[10];
  P1_min[1] = args[11];
  P2_plus[0] = args[12];
  P2_plus[1] = args[13];
  P2_min[0] = args[14];
  P2_min[1] = args[15];
}

#endif
