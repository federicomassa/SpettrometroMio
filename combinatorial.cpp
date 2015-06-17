#ifndef COMBINATORIAL_CPP
#define COMBINATORIAL_CPP

// combinatorial identification for two linear tracks in xz

class combinatorial {
 private:
  UInt_t fNdet; //number of detectors >= 3
  Double_t best_chi2;

  /* decimal representation of configuration:
     if = 0 ----> x1 = x1_init
     if = 2^(fNdet)-1 ---> x1 = x2_init
     else intermediate situations..see Config(UInt_t)
   */
  UInt_t best_config; 

  Double_t* z; //z array
  Double_t* x1_init;
  Double_t* x2_init;
  Double_t* x1; //best first x array
  Double_t* x2; //best second x array
  Double_t* y1_init;
  Double_t* y2_init;
  Double_t* y1; //best first y array
  Double_t* y2; //best second y array
  
  Double_t GetChi2(Double_t*, Double_t*);
  static void GetBinary(UInt_t, UInt_t*, UInt_t);
 public:
  combinatorial() {}
  void SetPoints(UInt_t, Double_t*, Double_t*, Double_t*, Double_t*, Double_t*); //z array, first x array, first y, second x array, second y in whatever order (1 <-> 2)
  void Config(UInt_t);
  void GetCombination(Double_t*, Double_t*);
  UInt_t GetBestConfig();
  Double_t GetBestChi2();

};

#endif
