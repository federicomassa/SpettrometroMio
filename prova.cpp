#include "String.h"
#include <iostream>
#include <fstream>
using namespace std;

void prova() {
  String str;
  string str2;
  ifstream in("prova.dat");
  getline(in, str);

  for (int i = 0; i < 5; i++)
    cout << str.GetElement(i) << endl;
}

