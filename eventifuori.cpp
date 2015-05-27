#include "event.h"
#include <iostream>

using namespace std;

void eventifuori() {

  event ev;
  ev.SetEvent(647, 647, 16.00000,100.58881,82.74116,0.00156,1.48798,17.84821,0.00723,4.62957, -0.00128, 0.01375, 25.00000, -0.00690, -0.06478, 25.00000, 0.00126, 0.03072, 35.00000, -0.01098, -0.13545, 35.00000);
  // ev.SetEvent(1,1, 9.83739, 100.49947, 32.16430, 0.00577, 5.95144, 68.33596, 0.00271, 2.80985, 0.08372, -0.02893, 25.00000, -0.03807, 0.01335, 25.00000, 0.13810, -0.04734, 35.00000, -0.06537, 0.02213, 35.00000);
  ev.kinematic_z_reco();
  cout << "Kin: " << ev.K_z_k << endl;
  cout << "CDA: " << ev.CDA_mean() << endl;

  Double_t final_coord[8];

  ev.IterationCycle(final_coord);
  cout << endl;
  cout << final_coord[0] << endl;
  cout << final_coord[1] << endl;
  cout << final_coord[2] << endl;
  cout << final_coord[3] << endl;
  cout << final_coord[4] << endl;
  cout << final_coord[5] << endl;
  cout << final_coord[6] << endl;
  cout << final_coord[7] << endl;

  Double_t pars[5] = {0,0,0,0,0};
  pars[0] = final_coord[0];
  pars[2] = final_coord[4];

  cout << "Z RECO: " << ev.z_k(pars) << endl;
}
