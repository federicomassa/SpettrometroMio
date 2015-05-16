#include <iostream>
#include <string>
#include <ctype.h>
#include <sstream>

void opt_split(const string option, string &splitted_opt, double &value) {

  string option_buff = "";
  string value_str;
  
  for (unsigned int i = 0; i < option.size() + 1; i++) {
    if (isalpha(option[i]) || option[i] == ' ') {
      option_buff += option[i];
    }
    else {
      value_str = option.substr(i+2,option.size() - (i+2));
      break;
    }
  }
  splitted_opt = option_buff.substr(0,option_buff.size()-1);
  stringstream(value_str) >> value;
}
