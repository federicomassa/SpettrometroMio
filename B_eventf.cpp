#ifndef B_EVENT_CPP
#define B_EVENT_CPP

#include "linef.h"
#include "helixf.h"
#include <TGraph2D.h>
#include <TCanvas.h>

class B_event {
private:
  line line_in, line_out;
  helix h;
  Float_t gamma, charge;
  Float_t p1[3], p2[3];
  Float_t theta_in, phi_in; //IN angles
  Float_t theta_out, phi_out; //OUT angles
  Float_t init_speed[3]; //initial speed vector
  Float_t time_interval; //time parameter value at the end of the magnetic field zone
  Float_t p3[3], p4[3];
  Float_t omega, radius;
  static Float_t L; 
  static Float_t Delta_z; //Distance between last detectors

public:
  B_event();
  B_event(Float_t*, Float_t*, Float_t, Float_t);
  void SetB_event(Float_t*, Float_t*, Float_t, Float_t);
  TGraph2D* GetEventGraph();
  Float_t GetTheta_in();
  Float_t GetTheta_out();
  Float_t GetPhi_in();
  Float_t GetPhi_out();
  Float_t GetTimeInterval();
  Float_t GetRadius();
};

#endif
