#ifndef HELIX_CPP
#define HELIX_CPP

class helix {
private:
  static Float_t mass;
  static Float_t B;
  
  Float_t init_point[3];
  Float_t init_speed[3];
  Float_t charge;
  Float_t gamma;
  Float_t omega;
  
public:
  helix();
  helix(Float_t*, Float_t*, Float_t, Float_t);
  void SetHelix(Float_t*, Float_t*, Float_t, Float_t);
  void GetInitCoord(Float_t*);
  void GetInitSpeed(Float_t*);
  Float_t GetCharge();
  Float_t GetGamma();
  Float_t GetOmega();
  Float_t GetRadius();
  Float_t GetParameterAtX(Float_t);
  void GetPoint(Float_t, Float_t*);
};
  
#endif
