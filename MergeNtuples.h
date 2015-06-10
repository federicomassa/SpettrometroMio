#ifndef MERGENTUPLES_H
#define MERGENTUPLES_H

#include <TNtupleD.h>
#include <iostream>
#include <string>
#include <sstream>

TNtupleD* MergeNtuples(TNtupleD* n1, TNtupleD* n2, const char* name = "merge", const char* title = "Merged Ntuple") {

  TNtupleD* n_merge;
  n_merge = 0;

  Int_t n_entries;
  
  if (n1->GetEntries() != n2->GetEntries()) {
    std::cout << "ERROR: different number of entries" << std::endl;
    return n_merge;
  }

  else n_entries = n1->GetEntries();

  Int_t n1_branches = n1->GetNbranches();
  Int_t n2_branches = n2->GetNbranches();

  string* branches_name1 = new string[n1_branches];
  string* branches_name2 = new string[n2_branches];

  stringstream ss;

  for (Int_t i = 0; i < n1_branches; i++) {
    ss << n1->GetListOfBranches()->At(i)->GetName();
    ss >> branches_name1[i];
    ss.clear();
      }

  for (Int_t i = 0; i < n2_branches; i++) {
    ss << n2->GetListOfBranches()->At(i)->GetName();
    ss >> branches_name2[i];
    ss.clear();
      }

  string branches_name3 = "";
  
  for (Int_t i = 0; i < n1_branches + n2_branches; i++) {
    if (i < n1_branches)
      branches_name3 = branches_name3 + branches_name1[i] + ":";
    else if (i != n1_branches + n2_branches - 1)
      branches_name3 = branches_name3 + branches_name2[i - n1_branches] + ":";
    else branches_name3 = branches_name3 + branches_name2[i - n1_branches];
  }

  n_merge = new TNtupleD(name, title, branches_name3.c_str());

  Double_t* arr1 = new Double_t[n1_branches];
  arr1 = (Double_t*) n1->GetArgs();

  Double_t* arr2 = new Double_t[n2_branches];
  arr2 = (Double_t*) n2->GetArgs();

  Double_t* arr3 = new Double_t[n1_branches + n2_branches];

  for (Int_t i = 0; i < n_entries; i++) {
    n1->GetEntry(i);
    n2->GetEntry(i);
    
    for (Int_t j = 0; j < n1_branches + n2_branches; j++) {
      arr3[j] = (j < n1_branches) ? arr1[j] : arr2[j - n1_branches];
    }
    
    n_merge->Fill(arr3);
    
  }

  return n_merge;

}

#endif
