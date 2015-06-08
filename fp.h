#ifndef FP_H
#define FP_H

#include "fp.cpp"

fp::fp() {}

void fp::SetFP(int (*f)(double)) {
  classfunc = f;
}

int fp::GetFValue(double val) {
  return classfunc(val);
}

#endif
