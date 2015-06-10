#include "MergeNtuples.h"

#include <TNtupleD.h>
#include <TFile.h>

void Merge() {

  TFile* in1 = new TFile("../Spettrometro_Files/reconstruction.root");
  TFile* in2 = new TFile("../Spettrometro_Files/p_reco.root");
  TFile* out = new TFile("../Spettrometro_Files/merge.root", "recreate");

  TNtupleD* n1 = (TNtupleD*) in1->Get("impr_meas");
  TNtupleD* n2 = (TNtupleD*) in2->Get("nt_out");

  TNtupleD* n3;
  n3 = MergeNtuples(n1,n2);
  
  n3->Write();
  out->Close();
}
