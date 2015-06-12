#ifndef ITERKINFITP_H
#define ITERKINFITP_H

#include "IterKinFitP.cpp"
#include <TMatrixD.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
 
IterKinFitP::IterKinFitP() {
  isInitialized = false;
  isFirstIteration = true;
  isFinal = true;

  fNvar = 0;
  fNconstr = 0;
  fNpar = 0;

  init_meas = 0;
  step_parameter = 0.5;
  threshold = 1E-5;
  MaxIterationNumber = 2000;
}

// Set constraint function, constraint derivative function, sigmas
void IterKinFitP::Initialize(UInt_t nvar, UInt_t nconstr, UInt_t npar, Double_t* init_m, Double_t* init_p, TMatrixD (*Phi_FCN)(Double_t*), TMatrixD (*B_FCN)(Double_t*), TMatrixD (*A_FCN)(Double_t*), Double_t* err) {

  fNvar = nvar;
  fNconstr = nconstr;
  fNpar = npar;

  if (fNpar >= fNconstr) std::cout << "ERROR: parameters >= constraints" << std::endl;

  init_meas = new Double_t[fNvar];
  init_pars = new Double_t[fNpar];

  final_meas = new Double_t[fNvar];
  final_pars = new Double_t[fNpar];

  for (UInt_t i = 0; i < fNvar; i++) {
    init_meas[i] = init_m[i];
    init_pars[i] = init_p[i];
  }
  
  isFirstIteration = true;
  isInitialized = true;

  Constr_FCN = Phi_FCN;
  Der_FCN = B_FCN;
  PDer_FCN = A_FCN;

  VarMatrix.ResizeTo(fNvar,fNvar);

  //initialize at zero
  for (UInt_t i = 0; i < fNvar; i++) {
    for (UInt_t j = 0; j < fNvar; j++) {
      VarMatrix(i,j) = 0.0;
    }
  }
  
  for (UInt_t i = 0; i < fNvar; i++) {
    VarMatrix(i,i) = err[i]*err[i];
  }

}

void IterKinFitP::Reset() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  isInitialized = false;
  isFirstIteration = true;

  fNvar = 0;
  fNconstr = 0;
  fNpar = 0;

  init_meas = 0;
  
  step_parameter = 0.5;
  threshold = 1E-5;
}

void IterKinFitP::SetStepParameter(Double_t step) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  step_parameter = step;
}

void IterKinFitP::SetThreshold(Double_t th) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  threshold = th;
}

void IterKinFitP::SetMaxIterationNumber(UInt_t max) {
  MaxIterationNumber = max;
}


void IterKinFitP::SetConstrVector(Double_t* var, Double_t* par) {
  ConstrVector = Constr_FCN(var,par);
}

TMatrixD IterKinFitP::GetVarianceMatrix() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return VarMatrix;
}

Double_t IterKinFitP::chi2(Double_t* var)  
 {
   
  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  Double_t chi = 0;
  
  for (UInt_t i = 0; i < fNvar; i++)
    chi += TMath::Power((var[i] - init_meas[i]),2)/VarMatrix(i,i);
  
  return chi;
  
 }

Double_t IterKinFitP::GetChi2(Double_t* var) {


  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
 
  return chi2(var);
  
}

void IterKinFitP::SetAMatrix() {
  AMatrix = PDer_FCN(curr_pars);
}

void IterKinFitP::SetBMatrix() {
  BMatrix = Der_FCN(curr_meas);
}


 // R = (D_t V D)
void IterKinFitP::SetRMatrix() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  

  TMatrixD BMatrix_t(fNconstr, fNvar);
  BMatrix_t.Transpose(BMatrix);

  TMatrixD Temp_Matrix(fNvar,fNconstr);
  Temp_Matrix.Mult(VarMatrix, BMatrix);

  TMatrixD R(fNconstr, fNconstr);
  R.Mult(BMatrix_t, Temp_Matrix);

  RMatrix = R;
  RinvMatrix = R.InvertFast();
 }

void IterKinFitP::SetPMatrix() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  TMatrixD P(fNpar, fNpar);

  TMatrixD Temp_Matrix1(fNconstr, fNpar);
  TMatrixD AMatrix_t(fNconstr, fNpar);
  AMatrix_t.Transpose(AMatrix);

  TMatrixD R_inv(fNconstr, fNconstr);
  
  if (RMatrix.Determinant() != 0) 
  R_inv = RMatrix.InvertFast();

  else return;

  Temp_Matrix1.Mult(R_inv, AMatrix_t);

  P.Mult(AMatrix,Temp_Matrix1);
    
  PMatrix = P;
  PinvMatrix = P.InvertFast();
}

    
 
Double_t IterKinFitP::GetRDeterminant() {
  return RMatrix.Determinant();
}

Double_t IterKinFitP::GetPDeterminant() {
  return PMatrix.Determinant();
}
   

TMatrixD IterKinFitP::SetDelta_a_Vector() {

  TMatrixD Delta_a_V(fNpar,1);
  TMatrixD Temp_Constr_V(fNconstr, 1);
  
  Temp_Constr_V.Mult(RinvMatrix, 

}

TMatrixD IterKinFitP::Delta_y_Vector (Double_t* var) {

 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   TMatrixD x_V(fNvar,1);

   TMatrixD BMatrix_t(fNconstr,fNvar);
   BMatrix_t.Transpose(BMatrix);

   if (RMatrix.Determinant() == 0) std::cout << "ERROR: R determinant = 0" << std::endl;

   TMatrixD Temp_Matrix1(fNvar,fNconstr);
   Temp_Matrix1.Mult(VarMatrix,BMatrix);
   
   Temp_Matrix1 *= R.InvertFast();
   
   TMatrixD Temp_Matrix2(fNconstr,1); //D_t x_G
   Temp_Matrix2.Mult(D_t,x_Guess(var));

   Temp_Matrix2 -= Constr_FCN(var);

   x_V.Mult(Temp_Matrix1,Temp_Matrix2);
   
   return x_V;
 }  

   
   
void IterKinFitP::Iterate(Double_t* old_var, Double_t* new_var) {
  
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

   TMatrixD x_G(fNvar,1);
   TMatrixD x_V(fNvar,1);
   TMatrixD D(fNvar,fNconstr);

   TMatrixD constr_V(fNconstr,1);

   x_G = x_Guess(old_var);
   constr_V = GetConstraintVector(old_var);

   x_V = x_Vector(old_var);
   for (UInt_t i = 0; i < fNvar; i++) {
     new_var[i] = old_var[i] + step_parameter*(x_V(i,0) - x_G(i,0));
   }


 }
   

void IterKinFitP::Minimize(Double_t* final_var) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   Double_t old_var[fNvar];
   Double_t new_var[fNvar];
   iteration_no = 0;     


     do {
       
       if (isFirstIteration) {
	 for (UInt_t i = 0; i < fNvar; i++) {
	   old_var[i] = init_meas[i];}
	 isFirstIteration = false;
       }
       
       Iterate(old_var, new_var);
       iteration_no++;

       //reset after check
       isFinal = true;
       
       //check if something's changed
       for (UInt_t i = 0; i < fNvar; i++) {
	 isFinal = isFinal && TMath::Abs((new_var[i] - old_var[i])/old_var[i]) < threshold;
       }
       
     
       for (UInt_t i = 0; i < fNvar; i++) old_var[i] = new_var[i];
     }
     
     while (!isFinal);
     
     iteration_no -= 1; //last iteration did not change anything above threshold
   
   for (UInt_t i = 0; i < fNvar; i++) {
     final_var[i] = new_var[i];
     final_meas[i] = new_var[i];
   }

 }


void IterKinFitP::Minimize(Double_t* final_var, TGraph* graph) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

 
   Double_t old_var[fNvar];
   Double_t new_var[fNvar];
   iteration_no = 0;     

   TGraph* g = new TGraph[fNvar];
   string graph_title = "Iteration plot of var ";
   string graph_name = "iter_graph";
   string graph_number;
   stringstream ss;

   for (UInt_t i = 0; i < fNvar; i++) {
     ss << i;
     ss >> graph_number;
     ss.clear();
     g[i].SetTitle((graph_title+graph_number).c_str());
     g[i].SetName((graph_name + graph_number).c_str());
   }


     do {
       
       if (isFirstIteration) {

	 for (UInt_t i = 0; i < fNvar; i++) {
	   g[i].SetPoint(0, 0.0, init_meas[i]);
	 }


	 for (UInt_t i = 0; i < fNvar; i++) {
	   old_var[i] = init_meas[i];}
	 
	 isFirstIteration = false;
       }
       
       Iterate(old_var, new_var);
       iteration_no++;


       //reset after check
       isFinal = true;
       
       //check if something's changed
       for (UInt_t i = 0; i < fNvar; i++) {
	 isFinal = isFinal && TMath::Abs((new_var[i] - old_var[i])/old_var[i]) < threshold;
       }
       
       if (!isFinal) {        // fill graph array
	 for (UInt_t i = 0; i < fNvar; i++) {
	   g[i].SetPoint(iteration_no, Double_t(iteration_no), new_var[i]);
	 }
       }
     
       for (UInt_t i = 0; i < fNvar; i++) old_var[i] = new_var[i];
     }
     
     while (!isFinal && iteration_no <= MaxIterationNumber);
     
     iteration_no -= 1; //last iteration did not change anything above threshold
   
   for (UInt_t i = 0; i < fNvar; i++) {
     final_var[i] = new_var[i];
     final_meas[i] = new_var[i];
     graph[i] = g[i];
   }

}

UInt_t IterKinFitP::GetIterationNumber() {
    if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

    return iteration_no;
}


void IterKinFitP::PrintResult() {

  std::cout << std::fixed << std::setprecision(10);

  std::cout << "---------------------------" << std::endl;
  std::cout << "Minimization result" << std::endl;
  std::cout << "---------------------------" << std::endl;

  std::cout << "Vars: " << fNvar << std::endl;
  std::cout << "Constr: " << fNconstr << std::endl;
  std::cout << "---------------------------" << std::endl;

  if (iteration_no != MaxIterationNumber)
  std::cout << "Converged after " << iteration_no << " iterations." << std::endl;
  else
    std::cout << "WARNING: Maximum iteration number reached." << std::endl;

  std::cout << "---------------------------" << std::endl;
  std::cout << "Chi2: " << '\t' << GetChi2(final_meas) << std::endl;
  std::cout << "---------------------------" << std::endl;

  for (UInt_t i = 0; i < fNvar; i++) {
    std::cout << "Var" << i << ": " << '\t' << final_meas[i] << std::endl;
  }

  std::cout << "---------------------------" << std::endl;
  
}

#endif 

    
