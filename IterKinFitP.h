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
  isRValid = true;
  isPValid = true;
  isValid = true;
  isOverMax = false;
  fNvar = 0;
  fNconstr = 0;
  fNpar = 0;

  init_meas = 0;
  step_parameter = 0.5;
  threshold = 1E-5;
  MaxIterationNumber = 2000;
}

// Set constraint function, constraint derivative function, sigmas
void IterKinFitP::Initialize(UInt_t nvar, UInt_t nconstr, UInt_t npar, Double_t* init_m, Double_t* init_p, TMatrixD (*Phi_FCN)(Double_t*, Double_t*), TMatrixD (*B_FCN)(Double_t*, Double_t*), TMatrixD (*A_FCN)(Double_t*, Double_t*), Double_t* err) {

  fNvar = nvar;
  fNconstr = nconstr;
  fNpar = npar;

  if (fNpar >= fNconstr) std::cout << "ERROR: parameters >= constraints" << std::endl;

  init_meas = new Double_t[fNvar];
  init_pars = new Double_t[fNpar];

  curr_meas = new Double_t[fNvar];
  curr_pars = new Double_t[fNpar];

  final_meas = new Double_t[fNvar];
  final_pars = new Double_t[fNpar];

  for (UInt_t i = 0; i < fNvar; i++) {
    init_meas[i] = init_m[i];
  }

  for (UInt_t i = 0; i < fNpar; i++) {
    init_pars[i] = init_p[i];
  }
  
  isFirstIteration = true;
  isInitialized = true;
  isRValid = true;
  isPValid = true;
  isValid = true;
  isOverMax = false;

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

Double_t IterKinFitP::chi2()  
 {
   
  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

  Double_t chi = 0;
  
  for (UInt_t i = 0; i < fNvar; i++)
    chi += TMath::Power((curr_meas[i] - init_meas[i]),2)/VarMatrix(i,i);
  
  return chi;
  
 }

Double_t IterKinFitP::GetChi2() {


  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
 
  return chi2();
  
}

void IterKinFitP::SetConstrVector() {
  ConstrVector.ResizeTo(fNconstr, 1);
  ConstrVector = Constr_FCN(curr_meas,curr_pars);
}

void IterKinFitP::SetAMatrix() {
  AMatrix.ResizeTo(fNpar, fNconstr);
  AMatrix = PDer_FCN(curr_meas, curr_pars);
}

void IterKinFitP::SetBMatrix() {

  BMatrix.ResizeTo(fNvar, fNconstr);
  BMatrix = Der_FCN(curr_meas, curr_pars);
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

  RMatrix.ResizeTo(fNconstr, fNconstr);
  RinvMatrix.ResizeTo(fNconstr, fNconstr);
  
  RMatrix = R;

  if (RMatrix.Determinant() != 0) 
    if (fNconstr > 1)
      RinvMatrix = RMatrix.InvertFast();
    else 
      RinvMatrix = RMatrix.Invert();
  
  else {isRValid = false; return;}
 }

void IterKinFitP::SetPMatrix() {

  if (isRValid) {
    
    if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
    
    TMatrixD P(fNpar, fNpar);
    
    TMatrixD Temp_Matrix1(fNconstr, fNpar);
    TMatrixD AMatrix_t(fNconstr, fNpar);
    AMatrix_t.Transpose(AMatrix);
    
    Temp_Matrix1.Mult(RinvMatrix, AMatrix_t);
    
    PMatrix.ResizeTo(fNpar, fNpar);
    PinvMatrix.ResizeTo(fNpar, fNpar);
    
    P.Mult(AMatrix,Temp_Matrix1);
    
    PMatrix = P;
    
    if (PMatrix.Determinant() != 0) 
      if (fNpar > 1)
	PinvMatrix = PMatrix.InvertFast();
      else 
	PinvMatrix = PMatrix.Invert();
    
    else {isPValid = false; return;}
  }

  else return;

  }
  
  
 
Double_t IterKinFitP::GetRDeterminant() {
  return RMatrix.Determinant();
}

Double_t IterKinFitP::GetPDeterminant() {
  return PMatrix.Determinant();
}
   

void IterKinFitP::SetDelta_a() {
  
  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

 Delta_a.ResizeTo(fNpar, 1);

 TMatrixD Temp_Matrix1(fNconstr, 1);
 TMatrixD Temp_Matrix2(fNpar, 1);
 TMatrixD Temp_Matrix3(fNpar, 1);

 Temp_Matrix1.Mult(RinvMatrix, ConstrVector);
 Temp_Matrix2.Mult(AMatrix, Temp_Matrix1);
 Temp_Matrix3.Mult(PinvMatrix, Temp_Matrix2);

 for (UInt_t i = 0; i < fNpar; i++)
   Temp_Matrix3(i,0) = -Temp_Matrix3(i,0);

 Delta_a = Temp_Matrix3;
  
}

void IterKinFitP::SetLambda() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  
  Lambda.ResizeTo(fNconstr, 1);
  
  TMatrixD Temp_Matrix(fNconstr,1);
  TMatrixD AMatrix_t(fNconstr, fNpar);
  
  AMatrix_t.Transpose(AMatrix);
  Temp_Matrix.Mult(AMatrix_t, Delta_a);
  
  TMatrixD Sum(fNconstr, 1);
  Sum.Plus(ConstrVector, Temp_Matrix);

  Lambda.Mult(RinvMatrix, Sum);
 }  


void IterKinFitP::SetDelta_y() {

  if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
  
  Delta_y.ResizeTo(fNvar, 1);
  
  TMatrixD TempVar(fNvar, fNvar);
  TMatrixD Temp_y(fNvar, 1);

  for (UInt_t i = 0; i < fNvar; i++) {
      for (UInt_t j = 0; j < fNvar; j++) {
	TempVar(i,j) = -VarMatrix(i,j);
      }
  }
  

  Temp_y.Mult(BMatrix, Lambda);
  Delta_y.Mult(TempVar, Temp_y);

 }  

void IterKinFitP::SetMatrices() {
  
       SetConstrVector();
       SetAMatrix();
       SetBMatrix();
       SetRMatrix();
       SetPMatrix();
       if (isRValid && isPValid) {
       SetDelta_a();
       SetLambda();
       SetDelta_y(); }
       else {isValid = false; return;}
}
   
   
void IterKinFitP::Iterate() {
  
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

 for (UInt_t i = 0; i < fNvar; i++) {
   curr_meas[i] += step_parameter*Delta_y(i,0);
 }

 for (UInt_t i = 0; i < fNpar; i++) {
   curr_pars[i] += step_parameter*Delta_a(i,0);
 }
 

}
   

void IterKinFitP::Minimize(Double_t* final_m, Double_t* final_p) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   Double_t old_var[fNvar];
   Double_t old_par[fNpar];

   iteration_no = 0;     

     do {

       
       if (isFirstIteration) {
	 for (UInt_t i = 0; i < fNvar; i++) {
	   curr_meas[i] = init_meas[i];
	   old_var[i] = curr_meas[i];
	 }

	 for (UInt_t i = 0; i < fNpar; i++) {
	   curr_pars[i] = init_pars[i];
	   old_par[i] = curr_pars[i];
	 }
	 
	 isFirstIteration = false;
       }
       
       // Computer matrices
       SetMatrices();
       if (isValid) {
	 Iterate();
	 iteration_no++;
	 
	 //reset after check
	 isFinal = true;
	 
	 //check if something's changed
	 for (UInt_t i = 0; i < fNvar; i++) {
	   isFinal = isFinal && TMath::Abs((curr_meas[i] - old_var[i])/old_var[i]) < threshold;
	 }
	 
	 for (UInt_t i = 0; i < fNpar; i++) {
	   isFinal = isFinal && TMath::Abs((curr_pars[i] - old_par[i])/old_par[i]) < threshold;
	 }
	 
	 for (UInt_t i = 0; i < fNvar; i++) old_var[i] = curr_meas[i];
	 for (UInt_t i = 0; i < fNpar; i++) old_par[i] = curr_pars[i];

       }
       else break;
       

     }
     
     while (!isFinal && iteration_no != MaxIterationNumber);
     
     if (iteration_no == MaxIterationNumber) isOverMax = true;
     
     iteration_no -= 1; //last iteration did not change anything above threshold
     
       for (UInt_t i = 0; i < fNvar; i++) {
	 final_meas[i] = curr_meas[i];
	 final_m[i] = curr_meas[i];
       }
       
       for (UInt_t i = 0; i < fNpar; i++) {
	 final_pars[i] = curr_pars[i];
	 final_p[i] = curr_pars[i];
       }
       
}

void IterKinFitP::Minimize(Double_t* final_m, Double_t* final_p, TGraph* graph) {
   
 if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;
   
   Double_t old_var[fNvar];
   Double_t old_par[fNpar];

   iteration_no = 0;     

   string graph_title_v = "Iteration plot of var ";
   string graph_name_v = "iter_graph_v";
   string graph_title_p = "Iteration plot of par ";
   string graph_name_p = "iter_graph_p";
   string graph_number;
   stringstream ss;

   for (UInt_t i = 0; i < fNvar; i++) {
     ss << i;
     ss >> graph_number;
     ss.clear();
     graph[i].SetTitle((graph_title_v+graph_number).c_str());
     graph[i].SetName((graph_name_v + graph_number).c_str());
   }

   for (UInt_t i = fNvar; i < fNvar + fNpar; i++) {
     ss << i;
     ss >> graph_number;
     ss.clear();
     graph[i].SetTitle((graph_title_p +graph_number).c_str());
     graph[i].SetName((graph_name_p + graph_number).c_str());
   }


     do {
       
       if (isFirstIteration) {
	 for (UInt_t i = 0; i < fNvar; i++) {
	   curr_meas[i] = init_meas[i];
	   old_var[i] = curr_meas[i];
	   graph[i].SetPoint(0, 0.0, curr_meas[i]);
	 }

	 for (UInt_t i = 0; i < fNpar; i++) {
	   curr_pars[i] = init_pars[i];
	   old_par[i] = curr_pars[i];
	   graph[i+fNvar].SetPoint(0, 0.0, curr_pars[i]);
	 }

	 isFirstIteration = false;
       }
       
       //Compute matrices
       SetMatrices();

       Iterate();
       iteration_no++;

       for (UInt_t i = 0; i < fNvar; i++) {
	 graph[i].SetPoint(iteration_no, Double_t(iteration_no), curr_meas[i]);
	 }

	 for (UInt_t i = 0; i < fNpar; i++) {
	   graph[i+fNvar].SetPoint(iteration_no, Double_t(iteration_no), curr_pars[i]);
	 }

       //reset after check
       isFinal = true;
       
       //check if something's changed
       for (UInt_t i = 0; i < fNvar; i++) {
	 isFinal = isFinal && TMath::Abs((curr_meas[i] - old_var[i])/old_var[i]) < threshold;
       } 
       
       for (UInt_t i = 0; i < fNpar; i++) {
	 isFinal = isFinal && TMath::Abs((curr_pars[i] - old_par[i])/old_par[i]) < threshold;
       }
       
       for (UInt_t i = 0; i < fNvar; i++) old_var[i] = curr_meas[i];
       for (UInt_t i = 0; i < fNpar; i++) old_par[i] = curr_pars[i];
     }
     
     while (!isFinal && iteration_no != MaxIterationNumber);

     if (iteration_no == MaxIterationNumber) isOverMax = true;

     iteration_no -= 1; //last iteration did not change anything above threshold
   
   for (UInt_t i = 0; i < fNvar; i++) {
     final_meas[i] = curr_meas[i];
     final_m[i] = curr_meas[i];
   }

   for (UInt_t i = 0; i < fNpar; i++) {
     final_pars[i] = curr_pars[i];
     final_p[i] = curr_pars[i];
   }

}



UInt_t IterKinFitP::GetIterationNumber() {
    if (!isInitialized) std::cout << "ERROR: First call Initialize" << std::endl;

    return iteration_no;
}

Bool_t IterKinFitP::IsOverMax() {
  return isOverMax;
}

void IterKinFitP::PrintResult() {

  std::cout << std::fixed << std::setprecision(10);

  std::cout << "---------------------------" << std::endl;
  std::cout << "Minimization result" << std::endl;
  std::cout << "---------------------------" << std::endl;

  std::cout << "Vars: " << fNvar << std::endl;
  std::cout << "Constr: " << fNconstr << std::endl;
  std::cout << "Pars: " << fNpar << std::endl;
  std::cout << "---------------------------" << std::endl;

  if (iteration_no != MaxIterationNumber)
  std::cout << "Converged after " << iteration_no << " iterations." << std::endl;
  else
    std::cout << "WARNING: Maximum iteration number reached." << std::endl;

  std::cout << "---------------------------" << std::endl;
  std::cout << "Chi2: " << '\t' << GetChi2(final_meas) << std::endl;
  std::cout << "---------------------------" << std::endl;

  for (UInt_t i = 0; i < fNvar; i++) {
    std::cout << "Var" << i << ": " << '\t' << init_meas[i] << '\t' << "->" << '\t' << final_meas[i] << std::endl;
  }

  for (UInt_t i = 0; i < fNpar; i++) {
    std::cout << "Par" << i << ": " << '\t' << init_pars[i] << '\t' << "->" << '\t' << final_pars[i] << std::endl;
  }

  std::cout << "---------------------------" << std::endl;
  
}

#endif 

    
