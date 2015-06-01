#ifndef B_EVENT_CPP
#define B_EVENT_CPP

#include "line.h"
#include "helix.h"
#include <TGraph2D.h>
#include <TCanvas.h>

class B_event {
private:
  line line_in, line_out;
  helix h;
  Double_t gamma, charge;
  Double_t p1[3], p2[3];
  Double_t theta_in, phi_in; //IN angles
  Double_t theta_out, phi_out; //OUT angles
  Double_t init_speed[3]; //initial speed vector
  Double_t time_interval; //time parameter value at the end of the magnetic field zone
  Double_t p3[3];
  Double_t omega, radius;
  static Double_t L; 

public:
  B_event();
  B_event(Double_t*, Double_t*, Double_t, Double_t);
  void SetB_event(Double_t*, Double_t*, Double_t, Double_t);
  TGraph2D* GetEventGraph();
  Double_t GetTheta_in();
  Double_t GetTheta_out();
  Double_t GetPhi_in();
  Double_t GetPhi_out();
  Double_t GetTimeInterval();
  Double_t GetRadius();
};

#endif
