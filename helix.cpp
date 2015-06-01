#ifndef HELIX_CPP
#define HELIX_CPP

class helix {
private:
  static Double_t mass;
  static Double_t B;
  
  Double_t init_point[3];
  Double_t init_speed[3];
  Double_t charge;
  Double_t gamma;
  Double_t omega;
  
public:
  helix();
  helix(Double_t*, Double_t*, Double_t, Double_t);
  void SetHelix(Double_t*, Double_t*, Double_t, Double_t);
  void GetInitCoord(Double_t*);
  void GetInitSpeed(Double_t*);
  Double_t GetCharge();
  Double_t GetGamma();
  Double_t GetOmega();
  Double_t GetRadius();
  Double_t GetParameter(Double_t*);
  void GetPoint(Double_t, Double_t*);
};
  
#endif
