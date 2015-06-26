#include <FGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <iostream>
void provaFGen() {

  FGenPhaseSpace event;
  Double_t mass[3] = {0.938, 0.139, 0.139};
  TLorentzVector beam(0,0,65,66);
  event.Initialize(beam, 3,mass);
  event.GenerateEvent();
  TLorentzVector* p = event.GetDecay(0);
  std::cout << p->Mag() << std::endl;
}
