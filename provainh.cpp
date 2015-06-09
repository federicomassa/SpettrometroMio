#include <iostream>

class Primitive {
private:
  int PRIV;
protected:
  int PROT;
public:
  int PUBL;
  void SET(int a, int b, int c) {PRIV = a; PROT = b; PUBL = c;}
  int PUBL_FCN() {return PRIV;}
};

class Derived : public Primitive {
private:
  int priv;
protected:
  int prot;
public:
  int publ;
  void set(int a, int b, int c) {priv = a; prot = b; publ = c;}
  int publ_fcn() {return PUBL;}
};

  void provainh() {

    Primitive* pr1 = new Primitive;
    pr1->SET(3,4,5);

    Derived* der1 = new Derived;
    der1->SET(6,7,8);
    
    std::cout << der1->publ_fcn() << std::endl;


    
    // std::cout << pr1->PUBL_FCN() << std::endl;
    // Derived* der1 = new Derived;
    // der1->set(6,7,8);
    // der1->SET(12,13,14);
    // der1->A.SET(0,0,0);
    
    // std::cout << der1->PUBL << std::endl;
    
    // std::cout << der1->publ << std::endl;

    // der1->PRIV = 1;
    // der1->priv = 15;

    // std::cout << der1->PRIV << std::endl;
    // std::cout << der1->publ_fcn() << std::endl;
    // std::cout << der1->priv << std::endl;

    // int* priv_p = &der1->priv;
    // *priv_p = 9;
    
    // std::cout << der1->priv << std::endl;

  }
