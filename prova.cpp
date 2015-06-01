#include "line.h"

#include <iostream>
#include <TMath.h>

using namespace std;

void prova() {

  double p0[3] = {0,0,0};
  double p1[3] = {1/TMath::Sqrt(3),1/TMath::Sqrt(3),1/TMath::Sqrt(3)};
  line line1;
  line1.SetLine(p0,p1);
  cout << "Theta: " << line1.GetTheta() << " Phi: " << line1.GetPhi() << endl;

  line line2;
  line2.SetLine(p1,p0);
  cout << "Theta: " << line2.GetTheta() << " Phi: " << line2.GetPhi() << endl;
  
  cout << "Check if equal: " << TMath::Pi() - line2.GetTheta() << '\t' << line2.GetPhi() - TMath::Pi() << endl;


}
  

