#include "fp.h"
#include "func.h"
#include <iostream>

void main_fp() {

  fp prova;
  prova.SetFP(func);
  std::cout <<  prova.GetFValue(7.8) << std::endl;

}
