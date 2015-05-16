#include <TRandom3.h>
#include <iostream>

void prova_seed() {
  TRandom3 rndgen;
  std::cout << rndgen.Uniform() << std::endl;
  UInt_t seed = rndgen.GetSeed();
  rndgen.SetSeed(seed);
}
