#ifndef ITERKINFIT_H
#define ITERKINFIT_H

#include "IterKinFit.cpp"
#include <TMatrixD.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <sstream>
 
IterKinFit::IterKinFit() {
  isDefined = false;
  isInitialized = false;
  isFirstIteration = true;
  isFinal = true;
  
  fNVar = 0;
  fNConstr = 0;
  init_meas = 0;
  step_parameter = 0.5;
  threshold = 1E-5;
}

void IterKinFit::SetIterKinFit(UInt_t nvar, UInt_t nconstr, Double_t* init) {
  fNVar = nvar;
  fNConstr = nconstr;
  init_meas = new Double_t[fNVar];
  final_meas = new Double_t[fNVar];
  for (UInt_t i = 0; i < fNVar; i++) {
    init_meas[i] = init[i];
  }
  
  isFirstIteration = true;
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
    VarMatrix(i,i) = err[i]*err[i];
  }

}

void IterKinFit::Reset() {

  isDefined = false;
  isInitialized = false;
  isFirstIteration = true;
  
  fNVar = 0;
  fNConstr = 0;
  init_meas = 0;
  
  step_parameter = 0.5;
}

void IterKinFit::SetStepParameter(Double_t step) {
  step_parameter = step;
}

void IterKinFit::SetThreshold(Double_t th) {
  threshold = th;
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
    chi += TMath::Power((var[i] - init_meas[i])/TMath::Sqrt(VarMatrix(i,i)),2);
  
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
  R.Mult(D_t, Temp_Matrix);

   return R;
 }

 // R = (D_t V D)
TMatrixD IterKinFit::GetRMatrix(TMatrixD D) {

  if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  

  TMatrixD D_t(fNConstr, fNVar);
  D_t.Transpose(D);

  TMatrixD Temp_Matrix(fNVar,fNConstr);
  Temp_Matrix.Mult(VarMatrix, D);

  TMatrixD R(fNConstr, fNConstr);
  R.Mult(D_t, Temp_Matrix);

   return R;
 }

Double_t IterKinFit::GetRDeterminant(Double_t* var) {

 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  TMatrixD R = GetRMatrix(var);

  return R.Determinant();
}

Double_t IterKinFit::GetRDeterminant(TMatrixD D) {

 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  TMatrixD R = GetRMatrix(D);

  return R.Determinant();
}



TMatrixD IterKinFit::x_Guess(Double_t* var) {

 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   TMatrixD x_G(fNVar,1);
   
   for (UInt_t i = 0; i < fNVar; i++) {
     x_G(i,0) = var[i]-init_meas[i];
   }

   return x_G;
}
   


TMatrixD IterKinFit::x_Vector (Double_t* var) {

 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   TMatrixD x_V(fNVar,1);
   TMatrixD D(fNVar, fNConstr);
   TMatrixD R(fNConstr,fNConstr);

   D = GetDerivativeMatrix(var);
   TMatrixD D_t(fNConstr,fNVar);
   D_t.Transpose(D);

   R = GetRMatrix(D);
   if (R.Determinant() == 0) std::cout << "ERROR: R determinant = 0" << std::endl;

   TMatrixD Temp_Matrix1(fNVar,fNConstr);
   Temp_Matrix1 = D;
   Temp_Matrix1 *= R.InvertFast();
   Temp_Matrix1.Mult(VarMatrix, Temp_Matrix1);

   TMatrixD Temp_Matrix2(fNConstr,1); //D_t x_G
   Temp_Matrix2.Mult(D_t,x_Guess(var));
   
   Temp_Matrix2 -= Constr_FCN(var);

   x_V.Mult(Temp_Matrix1,Temp_Matrix2);
   
   return x_V;
 }  

   
   
void IterKinFit::Iterate(Double_t* old_var, Double_t* new_var) {
  
 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

   TMatrixD x_G(fNVar,1);
   TMatrixD x_V(fNVar,1);
   TMatrixD D(fNVar,fNConstr);

   TMatrixD constr_V(fNConstr,1);

   x_G = x_Guess(old_var);
   constr_V = GetConstraintVector(old_var);

   x_V = x_Vector(old_var);
   for (UInt_t i = 0; i < fNVar; i++) {
     new_var[i] = old_var[i] + step_parameter*(x_V(i,0) - x_G(i,0));
   }


 }
   

void IterKinFit::Minimize(Double_t* final_var, UInt_t& iter_no) {
   
 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   Double_t old_var[fNVar];
   Double_t new_var[fNVar];
   iteration_no = 0;     


     do {
       
       if (isFirstIteration) {
	 for (UInt_t i = 0; i < fNVar; i++) {
	   old_var[i] = init_meas[i];}
	 isFirstIteration = false;
       }
       
       Iterate(old_var, new_var);
       iteration_no++;

       //reset after check
       isFinal = true;
       
       //check if something's changed
       for (UInt_t i = 0; i < fNVar; i++) {
	 isFinal = isFinal && TMath::Abs((new_var[i] - old_var[i])/old_var[i]) < threshold;
       }
       
     
       for (UInt_t i = 0; i < fNVar; i++) old_var[i] = new_var[i];
     }
     
     while (!isFinal);
     
     
   
   for (UInt_t i = 0; i < fNVar; i++) {
     final_var[i] = new_var[i];
   }

   iter_no = iteration_no;

 }


void IterKinFit::Minimize(Double_t* final_var, UInt_t& iter_no, TGraph* graph) {
   
 if (!isDefined) std::cout << "ERROR: First call SetIterKinFit" << std::endl;
  else if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   Double_t old_var[fNVar];
   Double_t new_var[fNVar];
   iteration_no = 0;     

   TGraph* g = new TGraph[fNVar];
   string graph_title;
   stringstream ss;

   for (UInt_t i = 0; i < fNVar; i++) {
     ss << i;
     ss >> graph_title;
     ss.clear();
     g[i].SetTitle(graph_title.c_str());
   }


     do {
       
       if (isFirstIteration) {

	 for (UInt_t i = 0; i < fNVar; i++) {
	   g[i].SetPoint(1, 0.0, init_meas[i]);
	 }


	 for (UInt_t i = 0; i < fNVar; i++) {
	   old_var[i] = init_meas[i];}
	 
	 isFirstIteration = false;
       }
       
       Iterate(old_var, new_var);
       iteration_no++;

       // fill graph array
       for (UInt_t i = 0; i < fNVar; i++) {
	 g[i].SetPoint(iteration_no+1, Double_t(iteration_no), new_var[i]);
       }

       //reset after check
       isFinal = true;
       
       //check if something's changed
       for (UInt_t i = 0; i < fNVar; i++) {
	 isFinal = isFinal && TMath::Abs((new_var[i] - old_var[i])/old_var[i]) < threshold;
       }
       
     
       for (UInt_t i = 0; i < fNVar; i++) old_var[i] = new_var[i];
     }
     
     while (!isFinal);
     
     
   
   for (UInt_t i = 0; i < fNVar; i++) {
     final_var[i] = new_var[i];
   }

   iter_no = iteration_no;

 }

#endif 

    
