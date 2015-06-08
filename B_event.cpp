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
  Double_t p1[3], p2[3]; //Before magnetic field points
  Double_t theta_in, phi_in; //IN angles
  Double_t theta_out, phi_out; //OUT angles
  Double_t init_speed[3]; //initial speed vector
  Double_t time_interval; //time parameter value at the end of the magnetic field zone
  Double_t p3[3], p4[3]; //After magnetic field points
  Double_t omega, radius, omegaT;

public:
  static Double_t L; 
  static Double_t Delta_z; //Distance between last detectors
  B_event();
  B_event(Double_t*, Double_t*, Double_t, Double_t);
  void SetB_event(Double_t*, Double_t*, Double_t, Double_t);
  TGraph2D* GetEventGraph();
  Double_t GetTheta_in();
  Double_t GetTheta_out();
  Double_t GetPhi_in();
  Double_t GetPhi_out();
  Double_t GetTimeInterval();
  Double_t GetOmega();
  Double_t GetOmegaT();
  Double_t GetRadius();
  void GetP3(Double_t*);
  void GetP4(Double_t*);
};

#endif
