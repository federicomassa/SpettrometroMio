#ifndef ITERKINFIT_H
#define ITERKINFIT_H

#include "IterKinFit.cpp"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <iostream>
 
IterKinFit::IterKinFit() {
  isDefined = false;
  isInitialized = false;
}

void IterKinFit::SetIterKinFit(Int_t nvar, Int_t nconstr, Double_t* init) {
  fNVar = nvar;
  fNConstr = nconstr;
  init_meas = new Double_t[fNVar];
  for (UInt_t i = 0; i < fNVar; i++) {
    init_meas[i] = init[i];
  }
  
  isDefined = true;
}

// Set constraint function, constraint derivative function, sigmas
void IterKinFit::Initialize(TMatrixD (*Phi_FCN)(Double_t*), TMatrixD (*D_FCN)(Double_t*), Double_t* err) {

  isInitialized = true;
  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;

  Constr_FCN = Phi_FCN;
  Der_FCN = D_FCN;

  TMatrixD VarMatrix(fNVar, fNVar);

  for (UInt_t i = 0; i < fNVar; i++) {
    VarMatrix(i,i) = err[i];
  }

}

TMatrixD IterKinFit::GetConstraintVector(Double_t* var) {
  
  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return Constr_FCN(var);

}

TMatrixD IterKinFit::GetDerivativeMatrix(Double_t* var) {

  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return Der_FCN(var);

}

TMatrixD IterKinFit::GetVarianceMatrix() {

  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return VarMatrix;
}

Double_t IterKinFit::chi2(Double_t* var)  
 {
   
  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  Double_t chi = 0;
  
  for (UInt_t i = 0; i < fNVar; i++)
    chi += TMath::Power((var[i] - init_meas[i])/sigma[i],2);
  
  return chi;
  
 }

Double_t IterKinFit::GetChi2(Double_t* var) {

  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
 
  return chi2(var);
  
}

 // R = (D_t V D)
TMatrixD IterKinFit::GetRMatrix(Double_t* var) {

  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  
  TMatrixD D(fNVar, fNConstr);
  D = Der_FCN(var);

  TMatrixD D_t(fNConstr, fNVar);
  D_t.Transpose(D);

  TMatrixD Temp_Matrix(fNVar,fNConstr);
  Temp_Matrix.Mult(VarMatrix, D);

  TMatrixD R(fNConstr, fNConstr);
  R.Mult(

  
   /* cout << "R: " << endl;  */
   /* cout << "    Det: " << R.Determinant() << endl; */
   
   /* R.Print(); */
 R.InvertFast();
   return R;
 }

 TMatrixD x_Guess(Double_t* measures_guess) {
   
   TMatrixD x_G(8,1);
   
   for (int i = 0; i < 8; i++) {
     x_G(i,0) = measures_guess[i]-measures[i];
   }

   return x_G;
 }
   


 TMatrixD x_Vector (TMatrixD D, TMatrixD x_G, TMatrixD Constr_Vector) {
   
   TMatrixD x_V(8,1);
   TMatrixD R(3,3);
   TMatrixD D_t(3,8);
   D_t.Transpose(D);

   R = GetRMatrix(D);

   TMatrixD A_1(8,3);
   A_1 = D;
   A_1 *= R;
   
   TMatrixD A_2(3,1); //D_t x_G
   A_2.Mult(D_t,x_G);
   
   A_2 -= Constr_Vector;

   x_V.Mult(A_1,A_2);
   
   return x_V;
 }  
   
 void Iterate(Double_t old_coord[], Double_t new_coord[]) {
  
   TMatrixD x_G(8,1);
   TMatrixD x_V(8,1);
   TMatrixD D(8,3);
   /* TMatrixD G_inv(8,8); */
   TMatrixD constr_V(3,1);


   if (isFirstIteration) {
     for (int i = 0; i < 8; i++) {
       old_coord[i] = measures[i];}
     isFirstIteration = false;
   }

   D = GetConstraintDerivativeMatrix(old_coord);
   /* cout << "D: " << endl; */
   /* D.Print(); */
   /* cout << "G_inv: " << endl; */
   /* G_inv = GetVarianceMatrix(); */
   /* G_inv.Print() */;
   
   x_G = x_Guess(old_coord);
   /* cout << "XG" << endl; */
   /* x_G.Print(); */
   /* cout << "FINE XG" << endl; */
   
   constr_V = GetConstraintsVector(old_coord);
   /* cout << "Constr_V: " << endl; */
   /* constr_V.Print(); */
   x_V = x_Vector(D, x_G, constr_V);
   /* cout << "x_V: " << endl; */
   /* x_V.Print(); */
   for (int i = 0; i < 8; i++) {
     new_coord[i] = old_coord[i] + 0.5*(x_V(i,0) - x_G(i,0));
   }


 }
   

 void IterationCycle(Double_t final_coord[]) {
   
   Double_t old_pars[5] = {0,0,0,0,0}, new_pars[5] = {0,0,0,0,0} /* final_pars[5] = {0,0,0,0,0} */ ;
   Double_t old_z, new_z;
   Double_t old_coord[8];
   Double_t new_coord[8];
   iteration_no = 0;
   /* ofstream fout("chi2_iter.dat"); */
   /* TFile* file_out = new TFile("iter.root", "recreate"); */
   /* TH1F* chi2_iter_hist = new TH1F("chi2_iter_hist", "Chi2 Iter Hist", 1000,0,1000); */
     /* chi2_iter_hist->SetBinContent(1,chi2_all(measures)); */

 

     for (int j = 0; j < 5; j++) par_constr[j] = 0;
     par_constr[0] = measures[0];
     par_constr[2] = measures[4];
     /* fout << chi2_all(measures) << '\t'; */
     /* fout << z_k(par_constr) << '\t'; */
     /* fout << measures[0] << '\t' << measures[1] << '\t' << measures[3] << endl; */
     /* chi2_iter_hist->SetBinContent(1,chi2_all(measures)); */
     

     do {
     Iterate(old_coord, new_coord);
     iteration_no++;

     old_pars[0] = old_coord[0];
     old_pars[2] = old_coord[4];

     new_pars[0] = new_coord[0];
     new_pars[2] = new_coord[4];

     old_z = z_k(old_pars);
     new_z = z_k(new_pars);
     
     //    cout << "It. no: " << iteration_no << "\t" << new_z << endl;

     for (int j = 0; j < 5; j++) par_constr[j] = 0;
     par_constr[0] = new_coord[0];
     par_constr[2] = new_coord[4];
     /* fout << chi2_all(new_coord) << '\t'; */
     /* fout << z_k(par_constr) << '\t'; */
     /* fout << new_coord[0] << '\t' << new_coord[1] << '\t' << new_coord[3] << endl; */
     /* chi2_iter_hist->SetBinContent(i+2,chi2_all(new_coord)); */

     
     for (int j = 0; j < 8; j++) old_coord[j] = new_coord[j];
     }
     
     while (TMath::Abs(new_z - old_z) > 1E-5 && iteration_no < 5000);
     
     //   cout << "-------------------------------------------" << endl;

     if (iteration_no == 5000) {/* cout << "MAX iteration number exceeded: decay number " << dec_no << endl;  */ isValid = false;
      }
     
   
   for (int i = 0; i < 8; i++) {
     final_coord[i] = new_coord[i];
   }

   /* final_pars[0] = final_coord[0]; */
   /* final_pars[2] = final_coord[4]; */
   
   /* if (iteration_no == 5000) cout << "z: " << z_k(final_pars) << endl; */


   /* chi2_iter_hist->Write(); */
   /* file_out->Close(); */
   /* fout.close(); */
 }

};


Double_t event::sigma = 0.001;

#endif 

    
