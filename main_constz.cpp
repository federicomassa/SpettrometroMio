#include <string>
#include <sstream>

#include "generation.cpp"
#include "detector.cpp"
#include "reconstruction.cpp"


void main_constz() {
  string gen_name, det_name, root_name;
  string gen_prefix = "events_z";
  string det_prefix = "measures_z";
  string root_prefix = "reconstruction_z";
  string str_buff;
  string str_dat = ".dat";
  string str_root = ".root";
  stringstream ss;
  double z_val;
  
  for (int i = 8; i < 10; i++) {
    z_val = double(i)*2;
    cout << "Starting analysis with z = " << z_val << "..." << endl;
    ss << z_val;
    ss >> str_buff;
    gen_name = gen_prefix + str_buff + str_dat;
    det_name = det_prefix + str_buff + str_dat;
    root_name = root_prefix + str_buff + str_root;
    
    ss.clear();

    cout << '\t' << "Generation..." << '\n' << endl;
    generation(gen_name.c_str(), "CONST Z");
    cout << '\t' << "Detector response..." << '\n' << endl;
    detector(gen_name.c_str(), det_name.c_str());
    cout << '\t' << "Reconstruction..." << '\n' << endl;
    reconstruction(gen_name.c_str(), det_name.c_str(), root_name.c_str());
  }


}
