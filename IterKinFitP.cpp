// Iterative Kinematic Fit with nonlinear constraints and unmeasured parameters
// for the version without parameters see IterKinFit class

// fNvar = number of measured variables
// fNconstr = number of constraints
// fNpar = number of unmeasured parameters

// B_ij = derivative of constraint j with respect to the i-th measured variable
// A_ij = derivative of constraint j with respect to the i-th unmeasured parameter
// V_ij = error matrix for measured variables

#ifndef ITERKINFITP_CPP
#define ITERKINFITP_CPP

#include <TMatrixD.h>
#include <TGraph.h>

class IterKinFitP { 
private:
  UInt_t fNvar, fNconstr, fNpar;
  Bool_t isInitialized, isFirstIteration, isFinal;
  Double_t step_parameter; //default = 0.5 (new_coord = old_coord + step_parameter*(x_vector-x_guess)
  Double_t threshold;
  UInt_t MaxIterationNumber;
  TMatrixD VarMatrix;
  TMatrixD ConstrVector;
  TMatrixD AMatrix; 
  TMatrixD BMatrix;
  TMatrixD RMatrix; //B_t V B
  TMatrixD RinvMatrix;
  TMatrixD PMatrix; // A R^-1 A_t
  TMatrixD PinvMatrix;
  TMatrixD Delta_y;
  TMatrixD Delta_a;
  TMatrixD Lambda;
  TMatrixD (*Constr_FCN)(Double_t*, Double_t*); //Constraint Vector Function
  TMatrixD (*Der_FCN)(Double_t*, Double_t*); //Derivative Matrix of measured vars
  TMatrixD (*PDer_FCN)(Double_t*, Double_t*); //Derivative Matrix of unmeasured pars

  Double_t* init_meas; //Measured
  Double_t* init_pars;

  Double_t* curr_meas; //While iterating
  Double_t* curr_pars;

  Double_t* final_meas; //Final measures
  Double_t* final_pars;

  UInt_t iteration_no;
  Double_t chi2();
  Double_t chi2(Double_t*);
  void SetConstrVector(); //to be set in order
  void SetAMatrix();
  void SetBMatrix();
  void SetRMatrix();
  void SetPMatrix();
  void SetDelta_a();
  void SetLambda();
  void SetDelta_y();
  void SetMatrices();
  void Iterate();
public:
  IterKinFitP();
  void Initialize(UInt_t, UInt_t, UInt_t, Double_t*, Double_t*, TMatrixD (*)(Double_t*, Double_t*), TMatrixD (*)(Double_t*, Double_t*), TMatrixD (*)(Double_t*, Double_t*), Double_t*);
  void Reset();
  void SetStepParameter(Double_t);
  void SetThreshold(Double_t);
  void SetMaxIterationNumber(UInt_t);
  TMatrixD GetVarianceMatrix();
  Double_t GetChi2();
  Double_t GetChi2(Double_t*);
  Double_t GetRDeterminant();
  Double_t GetPDeterminant();
  void Minimize(Double_t*, Double_t*);
  void Minimize(Double_t*, Double_t*, TGraph*);
  UInt_t GetIterationNumber();
  void PrintResult();

};

#endif
