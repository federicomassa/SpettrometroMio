#include "B_eventf.h"
#include <TGraph2D.h>
#include <iostream>

using namespace std;

void provaB() {

  B_event ev;
  Float_t p1[3] = {0,0,0};
  Float_t p2[3] = {0.001,0.002,1};
  Float_t gamma = 360;
  Float_t charge = 1;
  
  TGraph2D* g;

  ev.SetB_event(p1, p2, charge, gamma);
  g = ev.GetEventGraph();
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(7);
  g->Draw("P");

  cout << ev.GetTheta_in() << endl;
  cout << ev.GetPhi_in() << endl;
  cout << ev.GetTheta_in() - ev.GetTheta_out() << endl;
  cout << ev.GetPhi_in() - ev.GetPhi_out() << endl;
  cout << "TI: " << ev.GetTimeInterval() << endl;
  cout << "RAD: " << ev.GetRadius() << endl;
}
