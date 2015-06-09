#include <TFile.h>
#include <TGraph.h>

void ReadRootFile() {

  TFile* in = new TFile("provaIter.root"); //input file

  TGraph* gr = (TGraph*) in->Get("iter_graph0");
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerSize(2);
  gr->Draw("AP");

}
