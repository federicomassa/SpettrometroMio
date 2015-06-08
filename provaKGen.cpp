#include "KGen.h"
#include <iostream>

void provaKGen() {

  KGen ev;
  ev.Generate("CONST Z", 10.0);
  
  std::cout << ev.GetK_z() << std::endl;
  std::cout << ev.GetPi_plus_modp() << std::endl;
}
