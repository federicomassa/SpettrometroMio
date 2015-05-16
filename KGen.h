// Event generation class definition
// Cartesian system: z = K momentum direction

#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

class KGen {
 private:

  ////////// Physical and mathematical constants ////////////

  double c; // Light speed [m/ns] 
  double e; // Electron charge [aC]
  double pi; //3.14

  ///////////////// K parameters ////////////
  
  double K_tau0; // K mean lifetime [ns]
  double K_meanP; // K mean momentum [GeV]
  double K_sigmaP; // K ean momentum deviation [GeV]
  double K_Mass; // K mass [GeV]
  double pi_Mass; // Pi mass [GeV]

  /////////////// Generation  parameters ////////////
  
  double K_p; // K Momentum [GeV/c]
  double K_E; // K Energy [GeV]
  double K_lambda; // K Mean free path = gamma*beta*c*tau0 [m]
  double K_z; // K Decay point [m]
  double K_beta; 
  double K_gamma;
  
  double pi_plus_modp_star; // Positive Pi momentum in CMS [GeV/c]
  double pi_plus_E_star; // Positive Pi energy in CMS [GeV]
  double pi_min_modp_star; // Negative Pi momentum in CMS [GeV/c]
  double pi_min_E_star; // Negative Pi energy in CMS [GeV]
  double pi_plus_p_star[3]; //Positive Pi p vector in CMS [GeV/c]
  double pi_min_p_star[3]; //Negative Pi p vector in CMS [GeV/c]

  double pi_plus_modp; // Positive Pi momentum in LAB [GeV/c]
  double pi_plus_E; // Positive Pi energy in LAB [GeV]
  double pi_min_modp; // Negative Pi momentum in LAB [GeV/c]
  double pi_min_E; // Negative Pi energy in LAB [GeV]
  double pi_plus_p[3]; // Positive Pi p vector in LAB [GeV/c]
  double pi_min_p[3]; // Negative Pi p vector in LAB [GeV/c]

  double pi_plus_theta_star; // Positive Pi Polar theta in CMS [rad]
  double pi_min_theta_star; // Negative Pi Polar theta in CMS [rad]
  double pi_plus_phi_star; // Positive Pi Polar phi in CMS [rad]
  double pi_min_phi_star; // Positive Pi Polar phi in CMS [rad]

  double pi_plus_theta; // Positive Pi Polar theta in LAB [rad]
  double pi_min_theta; // Negative Pi Polar theta in LAB [rad]
  double pi_plus_phi; // Positive Pi Polar phi in LAB [rad]
  double pi_min_phi; // Positive Pi Polar phi in LAB [rad]


  /////////////// Random variables //////////////

  TRandom3 rndgen; //Random number generator.

 public:
  KGen() {
     ////////// Physical and mathematical constants ////////////

  c = 0.299792458; // Light speed [m/ns] 
  e = 0.160217657; // Electron charge [aC]
  pi = 4*TMath::ATan(1);

  ///////////////// K parameters ////////////
  
  K_tau0 = 0.08954; // K mean lifetime [ns]
  K_meanP = 100; // K mean momentum [GeV]
  K_sigmaP = 0.5; // K ean momentum deviation [GeV]
  K_Mass = 0.497614; // K mass [GeV]
  pi_Mass = 0.13957018; // Pi mass [GeV]
  
  K_z = -100;
  pi_plus_theta_star = -100;
  pi_min_theta_star = -100;
  pi_plus_phi_star = -100;
  pi_min_phi_star = -100;
  pi_plus_theta = -100;
  pi_min_theta = -100;
  pi_plus_phi = -100;
  pi_min_phi = -100;

  for (int i = 0; i < 3; i++) {
  pi_plus_p_star[i] = -100;
  pi_min_p_star[i] = -100;
  pi_plus_p[i] = -100;
  pi_min_p[i] = -100;
  }

  pi_plus_modp_star = -100;
  pi_min_modp_star = -100;
  pi_plus_modp = -100;
  pi_min_modp = -100;
  pi_plus_E_star = -100;
  pi_min_E_star = -100;
  pi_plus_E = -100;
  pi_min_E = -100;
  
  }

  void Generate(const char* gen_option = "", double option_val = 0)   // Event generation method
  {
    /*
      OPTIONS:
      CONST Z = <value>: generate events with fixed z
     */
    
    
    K_p = rndgen.Gaus(K_meanP, K_sigmaP);
    K_E = TMath::Sqrt(K_p*K_p + K_Mass*K_Mass);
    K_beta = 1 - 0.5*K_Mass*K_Mass/(K_p*K_p); // Second order Taylor expansion
    K_gamma = K_E/K_Mass;
    K_lambda = K_p/K_Mass*c*K_tau0;
    if (strcmp(gen_option,"") == 0)
      K_z = rndgen.Exp(K_lambda);
    else if(strcmp(gen_option,"CONST Z") == 0) {
      K_z = option_val;
    }

    
    
    pi_plus_theta_star = TMath::ACos(rndgen.Uniform(2)-1); //Cosine uniform distribution [-1,1]
    pi_plus_phi_star = rndgen.Uniform(2*pi); // Uniform distribution [0,2pi]
    pi_min_theta_star = pi - pi_plus_theta_star; //Opposite angle

    if (pi_plus_phi_star <= pi) pi_min_phi_star = pi_plus_phi_star + pi; //Opposite angle
    else pi_min_phi_star = pi_plus_phi_star - pi;

    pi_plus_phi = pi_plus_phi_star;
    pi_min_phi = pi_min_phi_star;

    pi_plus_E_star = K_Mass/2;
    pi_min_E_star = pi_plus_E_star;

    pi_plus_modp_star = TMath::Sqrt(pi_plus_E_star*pi_plus_E_star - pi_Mass*pi_Mass);
    pi_min_modp_star = pi_plus_modp_star;
    
    // CMS vectors fill
    pi_plus_p_star[0] = pi_plus_modp_star*TMath::Sin(pi_plus_theta_star)*TMath::Cos(pi_plus_phi_star); 
    pi_plus_p_star[1] = pi_plus_modp_star*TMath::Sin(pi_plus_theta_star)*TMath::Sin(pi_plus_phi_star);
    pi_plus_p_star[2] = pi_plus_modp_star*TMath::Cos(pi_plus_theta_star);

    pi_min_p_star[0] = -pi_plus_p_star[0];
    pi_min_p_star[1] = -pi_plus_p_star[1];
    pi_min_p_star[2] = -pi_plus_p_star[2];

    // LAB vectors fill (Boost)
    pi_plus_p[0] = pi_plus_p_star[0];
    pi_plus_p[1] = pi_plus_p_star[1];
    pi_plus_p[2] = K_gamma*(pi_plus_p_star[2] + K_beta*pi_plus_E_star);

    pi_min_p[0] = pi_min_p_star[0];
    pi_min_p[1] = pi_min_p_star[1];
    pi_min_p[2] = K_gamma*(pi_min_p_star[2] + K_beta*pi_min_E_star);
    


    pi_plus_modp = 0; // need to reset because the constructor is called just once in the main program
    pi_min_modp = 0;

    for (int i = 0; i < 3; i++) {
      pi_plus_modp += pi_plus_p[i]*pi_plus_p[i];
      pi_min_modp += pi_min_p[i]*pi_min_p[i];
  }
    
    pi_plus_E = TMath::Sqrt(pi_plus_modp + pi_Mass*pi_Mass); //p still squared
    pi_min_E = TMath::Sqrt(pi_min_modp + pi_Mass*pi_Mass);
    
    pi_plus_modp = TMath::Sqrt(pi_plus_modp); //now take the sqrt
    pi_min_modp = TMath::Sqrt(pi_min_modp);

    pi_plus_theta = TMath::ACos(pi_plus_p[2]/pi_plus_modp); // cos(theta) = p_z/|p|
    pi_min_theta = TMath::ACos(pi_min_p[2]/pi_min_modp);

  }

  void PrintEvent() { //Print event on screen
    cout << '\n' << endl;
    cout << "K generation..." << endl;

    cout << '\n' << endl;
    cout << '\t' << "Momentum:      " << K_p << endl;
    cout << '\t' << "Energy:        " << K_E << endl;
    cout << '\t' << "Momentum:      " << K_p << endl;
    cout << '\t' << "Beta:          " << K_beta << endl;
    cout << '\t' << "Gamma:         " << K_gamma << endl;
    cout << '\t' << "z of decay:    " << K_z << endl;
    cout << '\n' << endl;
    cout << "Pi generation..." << endl;
    cout << '\n' << endl;
    cout << '\t' << "CMS Pi+ E:     " << pi_plus_E_star << endl;
    cout << '\t' << "CMS Pi- E:     " << pi_min_E_star << endl;
    cout << '\t' << "CMS Pi+ p:     " << pi_plus_modp_star << endl;
    cout << '\t' << "CMS Pi- p:     " << pi_min_modp_star << endl;
    cout << '\t' << "CMS Pi+ px:    " << pi_plus_p_star[0] << endl;
    cout << '\t' << "CMS Pi+ py:    " << pi_plus_p_star[1] << endl;
    cout << '\t' << "CMS Pi+ pz:    " << pi_plus_p_star[2] << endl;
    cout << '\t' << "CMS Pi- px:    " << pi_min_p_star[0] << endl;
    cout << '\t' << "CMS Pi- py:    " << pi_min_p_star[1] << endl;
    cout << '\t' << "CMS Pi- pz:    " << pi_min_p_star[2] << endl;
    cout << '\t' << "CMS Pi+ theta: " << pi_plus_theta_star << endl;
    cout << '\t' << "CMS Pi- theta: " << pi_min_theta_star << endl;
    cout << '\t' << "CMS Pi+ phi:   " << pi_plus_phi_star << endl;
    cout << '\t' << "CMS Pi- phi:   " << pi_min_phi_star << endl;
    cout << '\t' << "LAB Pi+ E:     " << pi_plus_E << endl;
    cout << '\t' << "LAB Pi- E:     " << pi_min_E << endl;
    cout << '\t' << "LAB Pi+ p:     " << pi_plus_modp << endl;
    cout << '\t' << "LAB Pi- p:     " << pi_min_modp << endl;
    cout << '\t' << "LAB Pi+ px:    " << pi_plus_p[0] << endl;
    cout << '\t' << "LAB Pi+ py:    " << pi_plus_p[1] << endl;
    cout << '\t' << "LAB Pi+ pz:    " << pi_plus_p[2] << endl;
    cout << '\t' << "LAB Pi- px:    " << pi_min_p[0] << endl;
    cout << '\t' << "LAB Pi- py:    " << pi_min_p[1] << endl;
    cout << '\t' << "LAB Pi- pz:    " << pi_min_p[2] << endl;
    cout << '\t' << "LAB Pi+ theta: " << pi_plus_theta << endl;
    cout << '\t' << "LAB Pi- theta: " << pi_min_theta << endl;
    cout << '\t' << "LAB Pi+ phi:   " << pi_plus_phi << endl;
    cout << '\t' << "LAB Pi- phi:   " << pi_min_phi << endl;  
    cout << '\n' << endl;
  }
  
  // Writes event on a ASCII file. To view correctly, change your editor tab settings to at least 10 spaces
  void WriteEvent (ofstream &out, int ev_no = 0) { //out = output file, ev_no = event number
    if (ev_no == 1) out << "Event" << '\t'  << "z_dec" << '\t'  << "K_p" << '\t'  << "Pi+_p" << '\t'  << "Pi+_Th" << '\t'  << "Pi+_Phi" << '\t'  << "Pi-_p" << '\t'  << "Pi-_Th" << '\t'  << "Pi-_Phi" << '\n' << '\n'; 
    out << ev_no << fixed << setprecision(5) << '\t'  << K_z << '\t'  << K_p << '\t'  << pi_plus_modp << '\t'  << pi_plus_theta << '\t'  << pi_plus_phi << '\t'  << pi_min_modp << '\t'  << pi_min_theta << '\t'  << pi_min_phi << '\n'; 
  }

  //Functions to get event variables

  double GetK_z() { // Error if = -100
    double r;
    if (K_z != -100) r = K_z;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetK_p() { // Error if = -100
    double r;
    if (K_p != -100) r = K_p;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetK_E() { // Error if = -100
    double r;
    if (K_E != -100) r = K_E;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_plus_theta_star() { // Error if = -100
    double r;
    
    if (pi_plus_theta_star != -100) r = pi_plus_theta_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_min_theta_star() { // Error if = -100
    double r;
    
    if (pi_min_theta_star != -100) r = pi_min_theta_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_plus_phi_star() { // Error if = -100
    double r;
    
    if (pi_plus_phi_star != -100) r = pi_plus_phi_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_min_phi_star() { // Error if = -100
    double r;
    
    if (pi_min_phi_star != -100) r = pi_min_phi_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_plus_theta() { // Error if = -100
    double r;
    
    if (pi_plus_theta != -100) r = pi_plus_theta;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_min_theta() { // Error if = -100
    double r;
    
    if (pi_min_theta != -100) r = pi_min_theta;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_plus_phi() { // Error if = -100
    double r;
    
    if (pi_plus_phi != -100) r = pi_plus_phi;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_min_phi() { // Error if = -100
    double r;
    
    if (pi_min_phi != -100) r = pi_min_phi;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

  double GetPi_plus_px_star() {
    double r;
    if (pi_plus_p_star[0] != -100) r = pi_plus_p_star[0];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_plus_py_star() {
    double r;

    if (pi_plus_p_star[1] != -100) r = pi_plus_p_star[1];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_plus_pz_star() {
    double r;

    if (pi_plus_p_star[2] != -100) r = pi_plus_p_star[2];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};    
    
    return r;
  }

  double GetPi_min_px_star() {
    double r;
    if (pi_min_p_star[0] != -100) r = pi_min_p_star[0];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_min_py_star() {
    double r;

    if (pi_min_p_star[1] != -100) r = pi_min_p_star[1];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_min_pz_star() {
    double r;

    if (pi_min_p_star[2] != -100) r = pi_min_p_star[2];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};    
    
    return r;
  }

  double GetPi_plus_px() {
    double r;
    if (pi_plus_p[0] != -100) r = pi_plus_p[0];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_plus_py() {
    double r;

    if (pi_plus_p[1] != -100) r = pi_plus_p[1];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_plus_pz() {
    double r;

    if (pi_plus_p[2] != -100) r = pi_plus_p[2];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};    
    
    return r;
  }

  double GetPi_min_px() {
    double r;
    if (pi_min_p[0] != -100) r = pi_min_p[0];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_min_py() {
    double r;

    if (pi_min_p[1] != -100) r = pi_min_p[1];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};

    return r;
  }

  double GetPi_min_pz() {
    double r;

    if (pi_min_p[2] != -100) r = pi_min_p[2];
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};    
    
    return r;
  }


 double GetPi_plus_modp_star() { // Error if = -100
    double r;
    
    if (pi_plus_modp_star != -100) r = pi_plus_modp_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }
  
 double GetPi_min_modp_star() { // Error if = -100
    double r;
    
    if (pi_min_modp_star != -100) r = pi_min_modp_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

 double GetPi_plus_modp() { // Error if = -100
    double r;
    
    if (pi_plus_modp != -100) r = pi_plus_modp;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }
  
 double GetPi_min_modp() { // Error if = -100
    double r;
    
    if (pi_min_modp != -100) r = pi_min_modp;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }


 double GetPi_plus_E_star() { // Error if = -100
    double r;
    
    if (pi_plus_E_star != -100) r = pi_plus_E_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }
  
 double GetPi_min_E_star() { // Error if = -100
    double r;
    
    if (pi_min_E_star != -100) r = pi_min_E_star;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }

 double GetPi_plus_E() { // Error if = -100
    double r;
    
    if (pi_plus_E != -100) r = pi_plus_E;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }
  
 double GetPi_min_E() { // Error if = -100
    double r;
    
    if (pi_min_E != -100) r = pi_min_E;
    else {cout << "ERROR: Event not generated!" << endl; r = -100;};
    
    return r;
  }




};

