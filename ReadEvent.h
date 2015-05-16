#ifndef READEVENT_H
#define READEVENT_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;

void ReadEvent(ifstream &event_file, int &ev_no, double &K_z, double &K_p, double &pi_plus_modp, double &pi_plus_theta, double &pi_plus_phi, double &pi_min_modp, double &pi_min_theta, double &pi_min_phi) {
  string str;
  string entry = "";
  int j = 0;

  getline(event_file,str);

    for (unsigned int i = 0; i < str.size() + 1; i++) {
      
      if (!isalnum(str[i]) && str[i] != '.' && str[i] != '-') {
	j++;
	switch (j) {
	case(1):
	  stringstream(entry) >> ev_no;
	  entry = "";
	  break;
	case(2):
	  stringstream(entry) >> K_z;
	  entry = "";
	  break;
	case(3):
	  stringstream(entry) >> K_p;
	  entry = "";
	  break;
	case(4):
	  stringstream(entry) >> pi_plus_modp;
	  entry = "";
	  break;
	case(5):
	  stringstream(entry) >> pi_plus_theta;
	  entry = "";
	  break;
	case(6):
	  stringstream(entry) >> pi_plus_phi;
	  entry = "";
	  break;
	case(7):
	  stringstream(entry) >> pi_min_modp;
	  entry = "";
	  break;
	case(8):
	  stringstream(entry) >> pi_min_theta;
	  entry = "";
	  break;
	case(9):
	  stringstream(entry) >> pi_min_phi;
	  entry = "";
	  break;
	}
      }
      else {entry = entry + str.substr(i,1);}
      
      
      
   
    }
}


#endif
