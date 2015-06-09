#ifndef ITERKINFIT_H
#define ITERKINFIT_H

#include "IterKinFit.cpp"
#include <TMatrixD.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
 
IterKinFit::IterKinFit() {
  isInitialized = false;
  isFirstIteration = true;
  isFinal = true;
  
  fNVar = 0;
  fNConstr = 0;
  init_meas = 0;
  step_parameter = 0.2;
  threshold = 1E-5;
  MaxIterationNumber = 2000;
}

// Set constraint function, constraint derivative function, sigmas
void IterKinFit::Initialize(UInt_t nvar, UInt_t nconstr, Double_t* init, TMatrixD (*Phi_FCN)(Double_t*), TMatrixD (*D_FCN)(Double_t*), Double_t* err) {

  fNVar = nvar;
  fNConstr = nconstr;
  init_meas = new Double_t[fNVar];
  final_meas = new Double_t[fNVar];
  for (UInt_t i = 0; i < fNVar; i++) {
    init_meas[i] = init[i];
  }
  
  isFirstIteration = true;
  isInitialized = true;

  Constr_FCN = Phi_FCN;
  Der_FCN = D_FCN;

  VarMatrix.ResizeTo(fNVar,fNVar);

  //initialize at zero
  for (UInt_t i = 0; i < fNVar; i++) {
    for (UInt_t j = 0; j < fNVar; j++) {
      VarMatrix(i,j) = 0.0;
    }
  }
  
  for (UInt_t i = 0; i < fNVar; i++) {
    VarMatrix(i,i) = err[i]*err[i];
  }

}

void IterKinFit::Reset() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  isInitialized = false;
  isFirstIteration = true;
  
  fNVar = 0;
  fNConstr = 0;
  init_meas = 0;
  
  step_parameter = 0.5;
}

void IterKinFit::SetStepParameter(Double_t step) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  step_parameter = step;
}

void IterKinFit::SetThreshold(Double_t th) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  threshold = th;
}

void IterKinFit::SetMaxIterationNumber(UInt_t max) {
  MaxIterationNumber = max;
}

TMatrixD IterKinFit::GetConstraintVector(Double_t* var) {
  
  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return Constr_FCN(var);

}

TMatrixD IterKinFit::GetDerivativeMatrix(Double_t* var) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return Der_FCN(var);

}

TMatrixD IterKinFit::GetVarianceMatrix() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  return VarMatrix;
}

Double_t IterKinFit::chi2(Double_t* var)  
 {
   
  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  Double_t chi = 0;
  
  for (UInt_t i = 0; i < fNVar; i++)
    chi += TMath::Power((var[i] - init_meas[i]),2)/VarMatrix(i,i);
  
  return chi;
  
 }

Double_t IterKinFit::GetChi2(Double_t* var) {


  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
 
  return chi2(var);
  
}

 // R = (D_t V D)
TMatrixD IterKinFit::GetRMatrix(Double_t* var) {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  
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

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  

  TMatrixD D_t(fNConstr, fNVar);
  D_t.Transpose(D);
  
  TMatrixD Temp_Matrix(fNVar,fNConstr);
  Temp_Matrix.Mult(VarMatrix, D);

  TMatrixD R(fNConstr, fNConstr);
  R.Mult(D_t, Temp_Matrix);

   return R;
 }

Double_t IterKinFit::GetRDeterminant(Double_t* var) {

 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  TMatrixD R = GetRMatrix(var);

  return R.Determinant();
}

Double_t IterKinFit::GetRDeterminant(TMatrixD D) {

 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  TMatrixD R = GetRMatrix(D);

  return R.Determinant();
}



TMatrixD IterKinFit::x_Guess(Double_t* var) {

 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   TMatrixD x_G(fNVar,1);
   
   for (UInt_t i = 0; i < fNVar; i++) {
     x_G(i,0) = var[i]-init_meas[i];
   }

   return x_G;
}
   


TMatrixD IterKinFit::x_Vector (Double_t* var) {

 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   TMatrixD x_V(fNVar,1);
   TMatrixD D(fNVar, fNConstr);
   TMatrixD R(fNConstr,fNConstr);

   D = GetDerivativeMatrix(var);
   TMatrixD D_t(fNConstr,fNVar);
   D_t.Transpose(D);
   R = GetRMatrix(D);
   if (R.Determinant() == 0) std::cout << "ERROR: R determinant = 0" << std::endl;
   TMatrixD Temp_Matrix1(fNVar,fNConstr);
   Temp_Matrix1.Mult(VarMatrix,D);
   
   if (fNConstr > 1)
   Temp_Matrix1 *= R.InvertFast();
   else if(fNConstr == 1)
   Temp_Matrix1 *= R.Invert();
   
   TMatrixD Temp_Matrix2(fNConstr,1); //D_t x_G
   Temp_Matrix2.Mult(D_t,x_Guess(var));

   Temp_Matrix2 -= Constr_FCN(var);

   x_V.Mult(Temp_Matrix1,Temp_Matrix2);
   
   return x_V;
 }  

   
   
void IterKinFit::Iterate(Double_t* old_var, Double_t* new_var) {
  
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

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
   

void IterKinFit::Minimize(Double_t* final_var) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
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
     
     iteration_no -= 1; //last iteration did not change anything above threshold
   
   for (UInt_t i = 0; i < fNVar; i++) {
     final_var[i] = new_var[i];
     final_meas[i] = new_var[i];
   }

 }


void IterKinFit::Minimize(Double_t* final_var, TGraph* graph) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

 
   Double_t old_var[fNVar];
   Double_t new_var[fNVar];
   iteration_no = 0;     

   TGraph* g = new TGraph[fNVar];
   string graph_title = "Iteration plot of var ";
   string graph_name = "iter_graph";
   string graph_number;
   stringstream ss;

   for (UInt_t i = 0; i < fNVar; i++) {
     ss << i;
     ss >> graph_number;
     ss.clear();
     g[i].SetTitle((graph_title+graph_number).c_str());
     g[i].SetName((graph_name + graph_number).c_str());
   }


     do {
       
       if (isFirstIteration) {

	 for (UInt_t i = 0; i < fNVar; i++) {
	   g[i].SetPoint(0, 0.0, init_meas[i]);
	 }


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
       
       if (!isFinal) {        // fill graph array
	 for (UInt_t i = 0; i < fNVar; i++) {
	   g[i].SetPoint(iteration_no, Double_t(iteration_no), new_var[i]);
	 }
       }
     
       for (UInt_t i = 0; i < fNVar; i++) old_var[i] = new_var[i];
     }
     
     while (!isFinal && iteration_no <= MaxIterationNumber);
     
     iteration_no -= 1; //last iteration did not change anything above threshold
   
   for (UInt_t i = 0; i < fNVar; i++) {
     final_var[i] = new_var[i];
     final_meas[i] = new_var[i];
     graph[i] = g[i];
   }

}

UInt_t IterKinFit::GetIterationNumber() {
    if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

    return iteration_no;
}


void IterKinFit::PrintResult() {

  std::cout << std::fixed << std::setprecision(10);

  std::cout << "---------------------------" << std::endl;
  std::cout << "Minimization result" << std::endl;
  std::cout << "---------------------------" << std::endl;

  std::cout << "Vars: " << fNVar << std::endl;
  std::cout << "Constr: " << fNConstr << std::endl;
  std::cout << "---------------------------" << std::endl;

  if (iteration_no != MaxIterationNumber)
  std::cout << "Converged after " << iteration_no << " iterations." << std::endl;
  else
    std::cout << "WARNING: Maximum iteration number reached." << std::endl;

  std::cout << "---------------------------" << std::endl;
  std::cout << "Chi2: " << '\t' << GetChi2(final_meas) << std::endl;
  std::cout << "---------------------------" << std::endl;

  for (UInt_t i = 0; i < fNVar; i++) {
    std::cout << "Var" << i << ": " << '\t' << final_meas[i] << std::endl;
  }

  std::cout << "---------------------------" << std::endl;
  
}

#endif 

    
