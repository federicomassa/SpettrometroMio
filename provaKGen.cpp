#include "KGen.h"
#include <iostream>

void provaKGen() {

  KGen ev;
  ev.Generate("CONST Z", 10.0);
  

  std::cout << ev.GetPi_plus_px() << std::endl;
  std::cout << ev.GetPi_plus_py() << std::endl;
  std::cout << ev.GetPi_plus_pz() << std::endl;
  std::cout << ev.GetPi_min_px() << std::endl;
  std::cout << ev.GetPi_min_py() << std::endl;
  std::cout << ev.GetPi_min_pz() << std::endl;
}
