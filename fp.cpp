#ifndef FP_CPP
#define FP_CPP

class fp {
public:
  fp();
  int (*classfunc)(double);
  void SetFP(int (*)(double));
  int GetFValue(double);

};
  
#endif
