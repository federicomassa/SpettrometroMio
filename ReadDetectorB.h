#ifndef READDETECTORB_H
#define READDETECTORB_H

#include "GetElement.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;

void ReadDetectorB(ifstream &detector_file, Int_t &dec_no, Int_t &ev_no, Double_t &x1_plus, Double_t &y1_plus, Double_t &z1_plus, Double_t &x1_min, Double_t &y1_min, Double_t &z1_min, Double_t &x2_plus, Double_t &y2_plus, Double_t &z2_plus, Double_t &x2_min, Double_t &y2_min, Double_t &z2_min, Double_t &x3_plus, Double_t &y3_plus, Double_t &z3_plus, Double_t &x3_min, Double_t &y3_min, Double_t &z3_min, Double_t &x4_plus, Double_t &y4_plus, Double_t &z4_plus, Double_t &x4_min, Double_t &y4_min, Double_t &z4_min) {

  string str;
  

  getline(detector_file,str);

  dec_no = Int_t(GetElement(str,0));
  ev_no = Int_t(GetElement(str,1));
  x1_plus = GetElement(str,2);
  y1_plus = GetElement(str,3);
  z1_plus = GetElement(str,4);
  x1_min = GetElement(str,5);
  y1_min = GetElement(str,6);
  z1_min = GetElement(str,7);
  x2_plus = GetElement(str,8);
  y2_plus = GetElement(str,9);
  z2_plus = GetElement(str,10);
  x2_min = GetElement(str,11);
  y2_min = GetElement(str,12);
  z2_min = GetElement(str,13);
  x3_plus = GetElement(str,14);
  y3_plus = GetElement(str,15);
  z3_plus = GetElement(str,16);
  x3_min = GetElement(str,17);
  y3_min = GetElement(str,18);
  z3_min = GetElement(str,19);
  x4_plus = GetElement(str,20);
  y4_plus = GetElement(str,21);
  z4_plus = GetElement(str,22);
  x4_min = GetElement(str,23);
  y4_min = GetElement(str,24);
  z4_min = GetElement(str,25);

}

#endif
