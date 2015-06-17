#include "InitializeZPP.h"
#include "IterKinFitP.h"
#include "const_double.h"
#include <TFile.h>
#include <TNtupleD.h>

void iterkinfit_prova2() {

  TFile* in = new TFile("../Spettrometro_Files/measures.root");
  TNtupleD* nt = (TNtupleD*) in->Get("measures");
  Double_t* args = new Double_t[36];
  IterKinFitP* iter;
  Double_t* err = new Double_t[16];
  nt->GetEntry(0);
  args = nt->GetArgs();
  
  Double_t* init_measures = new Double_t[16];
  Double_t* init_pars = new Double_t[3];
  Double_t* final_measures = new Double_t[16];
  Double_t* final_pars = new Double_t[3];
  for (UInt_t i = 0; i < 16; i++) {
    init_measures[i] = args[20+i];
    err[i] = 0.001;
  }

  InitializeZPP(init_measures, init_pars);

  iter = new IterKinFitP;
  iter->Initialize(16,12,3,init_measures,init_pars, Constr, Der, PDer, err);
  iter->Minimize(final_measures, final_pars);
  iter->PrintResult();

}
