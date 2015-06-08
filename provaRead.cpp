#include <string>
#include <fstream>
#include <iostream>
#include "ReadDetectorB.h"

using namespace std;

void provaRead() {
  string str;
  ifstream in("../Spettrometro_Files/measuresB.dat");

  getline(in,str);
  getline(in,str);
  getline(in,str);
  getline(in,str);
  getline(in,str);
  getline(in,str);
  getline(in,str);

  Int_t dec_no, ev_no;
  Double_t
    x1_plus, y1_plus, z1_plus,
    x1_min, y1_min, z1_min,
    x2_plus, y2_plus, z2_plus,
    x2_min, y2_min, z2_min,
    x3_plus, y3_plus, z3_plus,
    x3_min, y3_min, z3_min,
    x4_plus, y4_plus, z4_plus,
    x4_min, y4_min, z4_min;

  ReadDetectorB(in, dec_no, ev_no, 
		x1_plus, y1_plus, z1_plus,
		x1_min, y1_min, z1_min,
		x2_plus, y2_plus, z2_plus,
		x2_min, y2_min, z2_min,
		x3_plus, y3_plus, z3_plus,
		x3_min, y3_min, z3_min,
		x4_plus, y4_plus, z4_plus,
		x4_min, y4_min, z4_min);

  cout << dec_no << endl;
  cout << ev_no << endl;
  cout << x1_plus << endl;
  cout << y1_plus << endl;
  cout << z1_plus << endl;
  cout << x1_min << endl;
  cout << y1_min << endl;
  cout << z1_min << endl;
  cout << x2_plus << endl;
  cout << y2_plus << endl;
  cout << z2_plus << endl;
  cout << x2_min << endl;
  cout << y2_min << endl;
  cout << z2_min << endl;
  cout << x3_plus << endl;
  cout << y3_plus << endl;
  cout << z3_plus << endl;
  cout << x3_min << endl;
  cout << y3_min << endl;
  cout << z3_min << endl;
  cout << x4_plus << endl;
  cout << y4_plus << endl;
  cout << z4_plus << endl;
  cout << x4_min << endl;
  cout << y4_min << endl;
  cout << z4_min << endl;

  
}
