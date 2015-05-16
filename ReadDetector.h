#ifndef READDETECTOR_H
#define READDETECTOR_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;

void ReadDetector(ifstream &detector_file, int &dec_no, int &ev_no, double &x1_plus, double &y1_plus, double &z1_plus, double &x1_min, double &y1_min, double &z1_min, double &x2_plus, double &y2_plus, double &z2_plus, double &x2_min, double &y2_min, double &z2_min) {
  string str;
  string entry = "";
  int j = 0;

  getline(detector_file,str);

    for (unsigned int i = 0; i < str.size() + 1; i++) {
      
      if (!isalnum(str[i]) && str[i] != '.' && str[i] != '-') {
	j++;
	switch (j) {
	case(1):
	  stringstream(entry) >> dec_no;
	  entry = "";
	  break;
	case(2):
	  stringstream(entry) >> ev_no;
	  entry = "";
	  break;
	case(3):
	  stringstream(entry) >> x1_plus;
	  entry = "";
	  break;
	case(4):
	  stringstream(entry) >> y1_plus;
	  entry = "";
	  break;
	case(5):
	  stringstream(entry) >> z1_plus;
	  entry = "";
	  break;
	case(6):
	  stringstream(entry) >> x1_min;
	  entry = "";
	  break;
	case(7):
	  stringstream(entry) >> y1_min;
	  entry = "";
	  break;
	case(8):
	  stringstream(entry) >> z1_min;
	  entry = "";
	  break;
	case(9):
	  stringstream(entry) >> x2_plus;
	  entry = "";
	  break;
	case(10):
	  stringstream(entry) >> y2_plus;
	  entry = "";
	  break;
	case(11):
	  stringstream(entry) >> z2_plus;
	  entry = "";
	  break;
	case(12):
	  stringstream(entry) >> x2_min;
	  entry = "";
	  break;
	case(13):
	  stringstream(entry) >> y2_min;
	  entry = "";
	  break;
	case(14):
	  stringstream(entry) >> z2_min;
	  entry = "";
	  break;
	}
      }
      else {entry = entry + str.substr(i,1);}
      
      
      
   
    }
}

#endif
