#ifndef ITERKINFIT_CPP
#define ITERKINFIT_CPP

#include <TMatrixD.h>
#include <TGraph.h>

class IterKinFit { 
private:
  UInt_t fNVar, fNConstr; //variables, constraint number
  Bool_t isDefined, isInitialized, isFirstIteration, isFinal;
  Double_t step_parameter; //default = 0.5 (new_coord = old_coord + step_parameter*(x_vector-x_guess)
  Double_t threshold;
  TMatrixD VarMatrix;
  TMatrixD (*Constr_FCN)(Double_t*); //Constraint Vector Function
  TMatrixD (*Der_FCN)(Double_t*); //Derivative Matrix function
  Double_t* init_meas;
  Double_t* final_meas;
  UInt_t iteration_no;
  Double_t chi2(Double_t*);
  TMatrixD GetRMatrix(Double_t*);
  TMatrixD GetRMatrix(TMatrixD);
  TMatrixD x_Guess(Double_t*);
  TMatrixD x_Vector(Double_t*);
  void Iterate(Double_t*, Double_t*);
public:
  IterKinFit();
  void SetIterKinFit(UInt_t, UInt_t, Double_t*);
  void Initialize(TMatrixD (*Phi_FCN)(Double_t*), TMatrixD (*D_FCN)(Double_t*), Double_t* err);
  void Reset();
  void SetStepParameter(Double_t);
  void SetThreshold(Double_t);
  TMatrixD GetConstraintVector(Double_t*);
  TMatrixD GetDerivativeMatrix(Double_t*);
  TMatrixD GetVarianceMatrix();
  Double_t GetChi2(Double_t*);
  Double_t GetRDeterminant(Double_t*);
  Double_t GetRDeterminant(TMatrixD);
  void Minimize(Double_t*, UInt_t&);
  void Minimize(Double_t*, UInt_t&, TGraph*);
  

};

#endif
