#ifndef GETELEMENT_H
#define GETELEMENT_H

#include <string>
#include <sstream>
#include <cstring>
#include <cctype>
#include <iostream>

using namespace std;

// n = 0 first element
Double_t GetElement(string str, unsigned int n) {
  
  if (str != "") {
  string str_buff = "";
  bool isValid = false;
  char ch_el[1];
  string str_el;
  stringstream ss;
  unsigned int j = 0;
  double return_val = 0;

  for (unsigned int i = 0; i < str.size() + 1; i++) {
    str_el = str.substr(i,1);

    // Previous check: list of allowed characters
    if (strcmp(str_el.c_str(), ".") != 0 && strcmp(str_el.c_str(), "-") != 0 && strcmp(str_el.c_str(), " ") != 0 && strcmp(str_el.c_str(), "\t") != 0 && !isdigit(str_el[0]) && strcmp(str_el.c_str(), "\0") != 0) continue;

    ss << str_el;
    ss >> ch_el;
   
    ss.clear();

    str_buff += str_el;
    
    // List of delimitators
    if (strcmp(str_el.c_str(), " ") == 0 || strcmp(str_el.c_str(), "\t") == 0 || strcmp(str_el.c_str(), "\0") == 0) {
      
	if (j == n) {isValid = true; stringstream(str_buff) >> return_val; break;}
	

	j++;
	str_buff = "";
    }
  }
 

  if (!isValid) cout << "ERROR: No matching element for String::GetElement(int)" << endl;

  isValid = false;  
  return return_val;}

  else return 0;
}

#endif
