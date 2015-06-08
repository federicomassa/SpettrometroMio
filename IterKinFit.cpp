#ifndef ITERKINFIT_CPP
#define ITERKINFIT_CPP

#include <TMatrixD.h>

class IterKinFit { 
private:
  UInt_t fNVar, fNConstr; //variables, constraint number
  Bool_t isDefined, isInitialized;
  TMatrixD VarMatrix;
  TMatrixD (*Constr_FCN)(Double_t*); //Constraint Vector Function
  TMatrixD (*Der_FCN)(Double_t*); //Derivative Matrix function
  Double_t* init_meas;
  Double_t chi2(Double_t*);
  TMatrixD GetRMatrix(Double_t*);
public:
  IterKinFit();
  void SetIterKinFit(Int_t, Int_t);
  void Initialize(TMatrixD (*Phi_FCN)(Double_t*), TMatrixD (*D_FCN)(Double_t*), Double_t* err);
  TMatrixD GetConstraintVector(Double_t*);
  TMatrixD GetDerivativeMatrix(Double_t*);
  TMatrixD GetVarianceMatrix();
  Double_t GetChi2(Double_t*);
  Double_t GetRDeterminant(Double_t*);
  
  
  

};

#endif
