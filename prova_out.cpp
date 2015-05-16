#include <iostream>
#include <fstream>

void prova_out() {

  ofstream out("prova.txt");
  double a = 3.2;
  double b = 4;
  out << a << b << '\t' << "ciao" << '\n' << endl;
  out.close();

}
