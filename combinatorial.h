#ifndef COMBINATORIAL_H
#define COMBINATORIAL_H

#include "IterKinFitP.h"
#include "const_double.h"
#include "InitializeZPP.h"
#include "combinatorial.cpp"
#include <TMath.h>
#include <iostream>

static Double_t sigma = 0.001;

void combinatorial::GetBinary(UInt_t decimal, UInt_t* binary, UInt_t nbits) {
  UInt_t dec;
  UInt_t bits; 

  dec = decimal;


  for (UInt_t i = 0; i < nbits; i++) {
    if (dec <= UInt_t(TMath::Nint(TMath::Power(2.0,Double_t(nbits-i))-1)) && dec >= UInt_t(TMath::Nint(TMath::Power(2.0,Double_t(nbits-1-i))))) binary[nbits-1 - i] = 1;
    else binary[nbits-1 - i] = 0;

    if (binary[nbits-1-i] == 1) dec -= TMath::Nint(TMath::Power(2.0,Double_t(nbits-1 - i)));
  }

}

Double_t combinatorial::GetChi2(Double_t* var1, Double_t* var2) { //equal sigmas

  Double_t chi = 0;

  for (UInt_t i = 0; i < fNdet; i++) { 
    chi += TMath::Power((var1[i]- var2[i])/sigma,2);
  }

  return chi;
}


void combinatorial::SetPoints(UInt_t det_number, Double_t* z_arr, Double_t* x_first, Double_t* y_first, Double_t* x_second, Double_t* y_second) {
  
  fNdet = det_number;
  z = new Double_t[fNdet];
  x1_init = new Double_t[fNdet];
  x2_init = new Double_t[fNdet];
  x1 = new Double_t[fNdet];
  x2 = new Double_t[fNdet];
  y1_init = new Double_t[fNdet];
  y2_init = new Double_t[fNdet];
  y1 = new Double_t[fNdet];
  y2 = new Double_t[fNdet];
  
  for (UInt_t i = 0; i < fNdet; i++) {
    z[i] = z_arr[i];
    x1_init[i] = x_first[i];
    y1_init[i] = y_first[i];
    x2_init[i] = x_second[i];
    y2_init[i] = y_second[i];
  }
}


// Given i as a decimal, converts it to binary and chooses the combination
void combinatorial::Config(UInt_t i) {


  UInt_t* binary = new UInt_t[fNdet];
  GetBinary(i,binary, fNdet); //config stores the analysed configuration
  
  for (UInt_t i = 0; i < fNdet; i++) {
    
    if (binary[i] == 0) {
      x1[i] = x1_init[i];
      y1[i] = y1_init[i];
      x2[i] = x2_init[i];
      y2[i] = y2_init[i];
    }
    else {
      x1[i] = x2_init[i]; 
      y1[i] = y2_init[i];
      x2[i] = x1_init[i];
      y2[i] = y1_init[i];
    }
   

  }

}

void combinatorial::GetCombination(Double_t* final_measures, Double_t* final_pars) {

  UInt_t* config = new UInt_t[fNdet];
  UInt_t best_config0 = 100;
  //improved measures array
  Double_t* x1_impr = new Double_t[fNdet];
  Double_t* x2_impr = new Double_t[fNdet];
  Double_t* y1_impr = new Double_t[fNdet];
  Double_t* y2_impr = new Double_t[fNdet];

  Double_t* init_m = new Double_t[16];
  Double_t* init_p = new Double_t[3];
  Double_t* final_m = new Double_t[16];
  Double_t* final_p = new Double_t[3];
  Double_t* err = new Double_t[16];

  for (UInt_t i = 0; i < 16; i++) err[i] = 0.001;

  Double_t chi2 = 100000000;
  Double_t chi2_buff = 10000000;
  IterKinFitP* iter;

  for (UInt_t i = 0; i < UInt_t(TMath::Nint((TMath::Power(2.0, Double_t(fNdet)))/2.0)); i++) {

  /* for (UInt_t i = 0; i < 1; i++) { */

    Config(i);
    
    for (UInt_t j = 0; j < 4; j++) 
      init_m[2*j] = x1[j];
    for (UInt_t j = 0; j < 4; j++) 
      init_m[2*j+1] = y1[j];
    for (UInt_t j = 0; j < 4; j++) 
      init_m[8 + 2*j] = x2[j];
    for (UInt_t j = 0; j < 4; j++) 
      init_m[8 + 2*j+1] = y2[j];

    InitializeZPP(init_m, init_p);
    

    //LS estimation of line passing through the fNdet points
    // x1 stores the first track, x2 stores the second
    /* { 
      Sx = Sz = Sxz = Szz = 0;
    for (UInt_t j = 0; j < fNdet; j++) {
      Sx += x1[j];
      Sz += z[j];
      Sxz += x1[j]*z[j];
      Szz += z[j]*z[j];
	
    }
   
    a[0] = (Double_t(fNdet)*Sxz - Sx*Sz)/(Double_t(fNdet)*Szz - Sz*Sz);
    b[0] = (Szz*Sx - Sxz*Sz)/(Double_t(fNdet)*Szz - Sz*Sz);

      Sx = Sz = Sxz = Szz = 0;
    for (UInt_t j = 0; j < fNdet; j++) {
      Sx += x2[j];
      Sz += z[j];
      Sxz += x2[j]*z[j];
      Szz += z[j]*z[j];
    }


    a[1] = (Double_t(fNdet)*Sxz - Sx*Sz)/(Double_t(fNdet)*Szz - Sz*Sz);
    b[1] = (Szz*Sx - Sxz*Sz)/(Double_t(fNdet)*Szz - Sz*Sz);
    

    for (UInt_t j = 0; j < fNdet; j++) {
      x1_impr[j] = a[0]*z[j] + b[0];
      x2_impr[j] = a[1]*z[j] + b[1];
    }
    
    chi2[0] = GetChi2(x1, x1_impr);
    chi2[1] = GetChi2(x2, x2_impr);

    chi2_mean_buff = chi2_mean;
    chi2_mean = (chi2[0] + chi2[1])/2.0;

    if (chi2_mean < chi2_mean_buff) {
      best_a1 = a[0];
      best_a2 = a[1];
      best_b1 = b[0];
      best_b2 = b[1];
      best_config0 = i;}
    else chi2_mean = chi2_mean_buff;
    
    } */

    ////////////////// ITERKINFIT //////////////////
    {
      iter = new IterKinFitP;
      iter->Initialize(16, 12, 3, init_m, init_p, Constr, Der, PDer, err);
      iter->SetThreshold(1E-7);
      iter->Minimize(final_m, final_p);

      chi2 = iter->GetChi2();

      /* std::cout << chi2_buff << '\t' << chi2 << std::endl; */
      /* std::cout << final_p[0] << '\t' << final_p[1] << '\t' << final_p[2] << std::endl; */
      /* std::cout << std::endl; */
      if (chi2 < chi2_buff) {
	best_config0 = i;
	chi2_buff = chi2;
	for (UInt_t k = 0; k < fNdet; k++) {
	  final_measures[k] = final_m[k];}
	for (UInt_t k = 0; k < 3; k++) {
	  final_pars[k] = final_p[k];}
      }
      
      delete iter;
    } //end of iterkinfit

   
    
  }

 /* std::cout << "\n////////////////////////\n"; */

  Config(best_config0); //x1, x2 now stores best configuration
  best_chi2 = chi2;
  best_config = best_config0;
  

}

Double_t combinatorial::GetBestChi2() {
  return best_chi2;
}

UInt_t combinatorial::GetBestConfig() {
  return best_config;
}  

/* void combinatorial::GetBest_ab(Double_t &a1, Double_t &b1, Double_t &a2, Double_t &b2) { */
  
/*   a1 = best_a1; */
/*   a2 = best_a2; */
/*   b1 = best_b1; */
/*   b2 = best_b2; */
  
/* } */

#endif
