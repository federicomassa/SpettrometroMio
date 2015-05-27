#include "reconstruction.cpp"

void reconstruction_all() {

  reconstruction("events_z0.dat", "measures_z0.dat", "reconstruction_z0.root");
  reconstruction("events_z2.dat", "measures_z2.dat", "reconstruction_z2.root");
  reconstruction("events_z4.dat", "measures_z4.dat", "reconstruction_z4.root");
  reconstruction("events_z6.dat", "measures_z6.dat", "reconstruction_z6.root");
  reconstruction("events_z8.dat", "measures_z8.dat", "reconstruction_z8.root");
  reconstruction("events_z10.dat", "measures_z10.dat", "reconstruction_z10.root");
  reconstruction("events_z12.dat", "measures_z12.dat", "reconstruction_z12.root");
  reconstruction("events_z14.dat", "measures_z14.dat", "reconstruction_z14.root");
  reconstruction("events_z16.dat", "measures_z16.dat", "reconstruction_z16.root");
  reconstruction("events_z18.dat", "measures_z18.dat", "reconstruction_z18.root");

}
