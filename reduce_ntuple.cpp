// takes ntuple with K_z K_p, p+, p-, theoretical points, measured points
// and converts it to ntuple without theoretical points

#include <TNtupleD.h>
#include <TFile.h>

void reduce_ntuple() {

  TFile* in = new TFile("../Spettrometro_Files/measures.root");
  TNtupleD* nt_old = (TNtupleD*) in->Get("measures");
  const char* branch_list = "K_z:K_p:p_plus:p_min:x1_plus:y1_plus:x2_plus:y2_plus:x3_plus:y3_plus:x4_plus:y4_plus:x1_min:y1_min:x2_min:y2_min:x3_min:y3_min:x4_min:y4_min";
  TNtupleD* nt_new = new TNtupleD("measures","Measures",branch_list);
  Double_t* args_old = new Double_t[36];
  Double_t* args_new = new Double_t[20];

  for (UInt_t i = 0; i < nt_old->GetEntries(); i++) {
    
  }

}
