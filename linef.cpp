#ifndef LINE_CPP
#define LINE_CPP

/*

PARAMETRIC 3D LINE

| x = init_point[0] + coeff[0]*t
| y = init_point[1] + coeff[1]*t
| z = init_point[2] + coeff[2]*t

coeff[0] = sin(theta)cos(phi)
coeff[1] = sin(theta)sin(phi)
coeff[2] = cos(theta)

 */

class line {
private:
  Float_t init_point[3]; 
  Float_t coeff[3];
  Float_t theta, phi;
public:
  line();
  line(Float_t*, Float_t*);
  line(Float_t*, Float_t, Float_t);
  void SetLine(Float_t*, Float_t*);
  void SetLine(Float_t*, Float_t, Float_t);
  Float_t GetParameterAtZ(Float_t);
  void GetPoint(Float_t, Float_t*);
  void GetInitCoords(Float_t*);
  void GetCoefficients(Float_t*);
  Float_t GetTheta();
  Float_t GetPhi();
};

#endif
