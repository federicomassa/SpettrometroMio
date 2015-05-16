#include <iostream>
#include <string>
#include <sstream>



void provastr() {
  const char c[] = "events_z15.dat";
  string str;
  stringstream ss;
  double val;

  ss << c;
  ss >> str;
  ss.clear();

  str = str.substr(8,2);
  ss << str;
  ss >> val;
  std::cout << val << std::endl;
}
