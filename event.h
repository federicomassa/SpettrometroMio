#ifndef EVENT_H
#define EVENT_H

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

class event {
 public:

static Double_t y2_plus_k_func(Double_t *par) {
  return par[2]/par[0]*par[1];
}

static Double_t x2_min_k_func(Double_t *par) {
  return par[3]/par[0]*par[2];
}

static Double_t y2_min_k_func(Double_t *par) {
  return par[2]/par[0]*par[4];
}

static Double_t z_k(Double_t *par) {
  return z1_plus - (z2_plus-z1_plus)/(par[2]-par[0])*par[0];
}

//x1+, y1+, x1-, y1-, analogous expression for the 2nd detector
 TMatrixD GetConstraintsVector(Double_t* coord) {
   TMatrixD constr_V(3,1);
   Double_t x1_p = coord[0];
   Double_t y1_p = coord[1];
   Double_t x1_m = coord[2];
   Double_t y1_m = coord[3];
   Double_t x2_p = coord[4];
   Double_t y2_p = coord[5];
   Double_t x2_m = coord[6];
   Double_t y2_m = coord[7];

   Double_t constr[3];
   constr[0] = y1_p*(x2_p - x1_p) - (y2_p - y1_p)*x1_p;
   constr[1] = y1_m*(x2_m - x1_m) - (y2_m - y1_m)*x1_m;
   constr[2] = x1_p*(x2_m - x1_m) - x1_m*(x2_p - x1_p);

   constr_V(0,0) = constr[0];
   constr_V(1,0) = constr[1];
   constr_V(2,0) = constr[2];

   return constr_V;

 }

  static Int_t dec_no, ev_no;
  Int_t iteration_no;
  static Double_t K_z, K_p, 
    pi_plus_modp, pi_plus_theta, pi_plus_phi,
    pi_min_modp, pi_min_theta, pi_min_phi,
    x1_plus, y1_plus, z1_plus, x1_min, y1_min, z1_min,
    x2_plus, y2_plus, z2_plus, x2_min, y2_min, z2_min,
    x1_plus_k, y1_plus_k, z1_plus_k, x1_min_k, y1_min_k, z1_min_k,
    x2_plus_k, y2_plus_k, z2_plus_k, x2_min_k, y2_min_k, z2_min_k,
    K_z_k, sigma;
  Double_t measures[8];
  Double_t coord_init[8];
  Bool_t isFirstIteration;
  Bool_t isValid;

  event() {
    dec_no = ev_no =
    K_z =  K_p =  
    pi_plus_modp =  pi_plus_theta =  pi_plus_phi = 
    pi_min_modp =  pi_min_theta =  pi_min_phi = 
    x1_plus =  y1_plus =  z1_plus =  x1_min =  y1_min =  z1_min = 
    x2_plus =  y2_plus =  z2_plus =  x2_min =  y2_min =  z2_min =
    x1_plus_k = y1_plus_k = z1_plus_k = x1_min_k = y1_min_k = z1_min_k =
    x2_plus_k = y2_plus_k = z2_plus_k = x2_min_k = y2_min_k = z2_min_k =
      K_z_k = -100;
    isFirstIteration = true;
  }
    
  void SetEvent(Int_t dec, Int_t ev, Double_t Kz, Double_t Kp, 
       Double_t piplusmodp, Double_t piplustheta, Double_t piplusphi, 
       Double_t piminmodp, Double_t pimintheta, Double_t piminphi, 
       Double_t x1plus, Double_t y1plus, Double_t z1plus, Double_t x1min, Double_t y1min, Double_t z1min,
		Double_t x2plus, Double_t y2plus, Double_t z2plus, Double_t x2min, Double_t y2min, Double_t z2min) {
    
    if (Kz > z1plus || Kz > z1min || Kz > z2plus || Kz > z2min)
      cout << "ERROR: decay after detector" << endl;

    dec_no = dec;
    ev_no = ev;
    K_z = Kz;
    K_p = Kp;
    pi_plus_modp = piplusmodp;
    pi_plus_theta = piplustheta;
    pi_plus_phi = piplusphi;
    pi_min_modp = piminmodp;
    pi_min_theta = pimintheta;
    pi_min_phi = piminphi;
    
    x1_plus = x1plus;
    measures[0] = x1plus;
    coord_init[0] = (z1plus - Kz)*TMath::Tan(piplustheta)*TMath::Cos(piplusphi);

    y1_plus = y1plus;
    measures[1] = y1plus;
    coord_init[1] = (z1plus - Kz)*TMath::Tan(piplustheta)*TMath::Sin(piplusphi);

    z1_plus = z1plus;

    x1_min = x1min;
    measures[2] = x1min;
    coord_init[2] = (z1min - Kz)*TMath::Tan(pimintheta)*TMath::Cos(piminphi);

    y1_min = y1min;
    measures[3] = y1min;
    coord_init[3] = (z1min - Kz)*TMath::Tan(pimintheta)*TMath::Sin(piminphi);


    z1_min = z1min;

    x2_plus = x2plus;
    measures[4] = x2plus;
    coord_init[4] = (z2plus - Kz)*TMath::Tan(piplustheta)*TMath::Cos(piplusphi);

    y2_plus = y2plus;
    measures[5] = y2plus;
    coord_init[5] = (z2plus - Kz)*TMath::Tan(piplustheta)*TMath::Sin(piplusphi);

    z2_plus = z2plus;

    x2_min = x2min;
    measures[6] = x2min;
    coord_init[6] = (z2min - Kz)*TMath::Tan(pimintheta)*TMath::Cos(piminphi);

    y2_min = y2min;
    measures[7] = y2min;
    coord_init[7] = (z2min - Kz)*TMath::Tan(pimintheta)*TMath::Sin(piminphi);

    z2_min = z2min;
    
    isFirstIteration = true;
    isValid = true;

  }



  
    

  void Reconstruction(TH1F* CDA_plus_hist, TH1F* CDA_min_hist, TH1F* CDA_mean_hist, TH1F* kinematic_z_hist, TH1F* iteration_hist, TH1F* iteration_cut_hist, TH1F* iteration_no_hist, TH2F* z_itno_iter_corr_hist, TH1F* kin_iter_diff_hist, TH2F* kin_iter_zz_hist, TH1F* chi2_init_hist, TH2F* chi2_corr_hist, TH1F* chi2_kin_hist, TH1F* chi2_iter_hist, TH2F* z_theta_corr_hist) {
    //Types of reconstruction: 1) Using the CDA of just one of the particles (pi+ or pi-) 2) Using both CDA
    Double_t chi2_val;
    Int_t npar = 5;
    Double_t* gin = new Double_t;
    Int_t iflag = 0;

    chi2_init_hist->Fill(chi2_all(coord_init));
    CDA_plus_hist->Fill(CDA_plus() - K_z);
    CDA_min_hist->Fill(CDA_min() - K_z);
    CDA_mean_hist->Fill(CDA_mean() - K_z);
    
    kinematic_z_reco();

    kinematic_z_hist->Fill(K_z_k - K_z);

    Double_t final_coord[8], pars[5];
    IterationCycle(final_coord);
    
    pars[0] = final_coord[0];
    pars[1] = final_coord[1];
    pars[2] = final_coord[4];
    pars[3] = final_coord[2];
    pars[4] = final_coord[3];
    
    //iteration number cut and chi2 cut (3 ndf, 0.1% probability)
    if (isValid) iteration_hist->Fill(z_k(pars) - K_z); 
    if (isValid && chi2_all(final_coord) <= 16.2662) iteration_cut_hist->Fill(z_k(pars) - K_z); 
    chi2_corr_hist->Fill(K_z - z_k(pars), chi2_all(final_coord));
    if (TMath::Abs(z_k(pars) - K_z) > 5) chi2_iter_hist->Fill(chi2_all(final_coord));
    chi2(npar, gin, chi2_val, pars, iflag);
    chi2_kin_hist->Fill(chi2_val);
    iteration_no_hist->Fill(iteration_no);

    kin_iter_diff_hist->Fill(K_z_k - z_k(pars));
    kin_iter_zz_hist->Fill(K_z - K_z_k, K_z - z_k(pars)); 
    if (TMath::Abs(K_z - K_z_k) < 0.01 && TMath::Abs(K_z - z_k(pars)) > 5)
      cout << "WARNING: Good Kin vs. Bad Iter at decay number: " << dec_no << endl;

    z_theta_corr_hist->Fill(K_z - z_k(pars), pi_plus_theta < pi_min_theta ? pi_plus_theta : pi_min_theta); 
    z_itno_iter_corr_hist->Fill(K_z - z_k(pars), double(iteration_no));
  }

  Double_t CDA_plus() {
    Double_t r;
    r = z1_plus - (z2_plus - z1_plus)*(x1_plus*(x2_plus-x1_plus) + y1_plus*(y2_plus-y1_plus))/((x2_plus - x1_plus)*(x2_plus - x1_plus) + (y2_plus - y1_plus)*(y2_plus - y1_plus));
    
    return r;
  }
   
  Double_t CDA_min() {
    Double_t r;
    r = z1_min - (z2_min - z1_min)*(x1_min*(x2_min-x1_min) + y1_min*(y2_min-y1_min))/((x2_min - x1_min)*(x2_min - x1_min) + (y2_min - y1_min)*(y2_min - y1_min));
    
    return r;
  }

  Double_t CDA_mean() {
    Double_t r;
    r = 0.5*(CDA_plus() + CDA_min());

    return r;
  }

static void chi2(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag)  
 {
  
  ///////// Parameters: x1+, y1+, x2+, x1-, y1-. The rest has been removed due to constraints.
  Double_t chi;
  

  chi = 
    TMath::Power((par[0] - x1_plus)/sigma,2) + 
    TMath::Power((par[1] - y1_plus)/sigma,2) + 
    TMath::Power((par[2] - x2_plus)/sigma,2) + 
    TMath::Power((y2_plus_k_func(par) - y2_plus)/sigma,2) +
    TMath::Power((par[3] - x1_min)/sigma,2) + 
    TMath::Power((par[4] - y1_min)/sigma,2) + 
    TMath::Power((x2_min_k_func(par) - x2_min)/sigma,2) + 
    TMath::Power((y2_min_k_func(par) - y2_min)/sigma,2); 
  
  f = chi + npar*0 + Double_t(iflag)*0 + (*gin)*0; //npar*0 to cancel compilation warnings
  
}

 Double_t chi2_all(Double_t coord[])  
 {
  
  
  Double_t chi;
  

  chi = 
    TMath::Power((coord[0] - x1_plus)/sigma,2) + 
    TMath::Power((coord[1] - y1_plus)/sigma,2) + 
    TMath::Power((coord[2] - x1_min)/sigma,2) + 
    TMath::Power((coord[3] - y1_min)/sigma,2) +
    TMath::Power((coord[4] - x2_plus)/sigma,2) + 
    TMath::Power((coord[5] - y2_plus)/sigma,2) + 
    TMath::Power((coord[6] - x2_min)/sigma,2) + 
    TMath::Power((coord[7] - y2_min)/sigma,2); 
 
  return chi;
  
}



void kinematic_z_reco() {

  Double_t final_pars[5], final_pars_err[5]; //Minimization results

  TMinuit* min = new TMinuit(5);
  min->SetPrintLevel(-1); //quiet mode
  min->SetFCN(chi2); //Set function to minimize

  Double_t arglist[10]; // ?
  Int_t ierflg = 0; // 0 if minimization converged succesfully

  arglist[0] = 1; // ?
  min->mnexcm("SET ERR", arglist, 1, ierflg); //Set err to 1 (standard value for chi2 minimization)

  // Initialize
  Double_t vstart[5];
  vstart[0] = x1_plus;
  vstart[1] = y1_plus;
  vstart[2] = x2_plus;
  vstart[3] = x1_min;
  vstart[4] = y1_min;
  
  //Initial parameter steps
  Double_t step[5] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001}; 
  
  // Define parameters
  min->mnparm(0, "x1+", vstart[0], step[0], 0,0,ierflg);
  min->mnparm(1, "y1+", vstart[1], step[1], 0,0,ierflg);
  min->mnparm(2, "x2+", vstart[2], step[2], 0,0,ierflg);
  min->mnparm(3, "x1-", vstart[3], step[3], 0,0,ierflg);
  min->mnparm(4, "y1-", vstart[4], step[4], 0,0,ierflg);

  arglist[0] = 500;
  arglist[1] = 1;
  
  // Minimization
  min->mnexcm("MIGRAD", arglist, 2, ierflg); 

   /* Double_t amin, edm, errdef;  */
   /* Int_t nvpar, nparx, icstat;  */

   /* // Show results  */
   /* min->mnstat(amin,edm,errdef,nvpar,nparx,icstat);  */

 
  
  min->GetParameter(0, final_pars[0], final_pars_err[0]); 
  min->GetParameter(1, final_pars[1], final_pars_err[1]);
  min->GetParameter(2, final_pars[2], final_pars_err[2]);
  min->GetParameter(3, final_pars[3], final_pars_err[3]);
  min->GetParameter(4, final_pars[4], final_pars_err[4]);

  x1_plus_k = final_pars[0];
  y1_plus_k = final_pars[1];
  x2_plus_k = final_pars[2];
  y2_plus_k = y2_plus_k_func(final_pars);
  
  x1_min_k = final_pars[3];
  y1_min_k = final_pars[4];
  x2_min_k = x2_min_k_func(final_pars);
  y2_min_k = y2_min_k_func(final_pars);

  K_z_k = z_k(final_pars);

  /* if (K_z_k - K_z > 4) cout << "Dec. no: " << dec_no << '\n' << "Ev. no: " << ev_no << endl; */

  /* cout << "FINAL RESULTS " << '\n'  */
  /*   << "x1+" << '\t' << P1_plus_reco[0] << '\n' */
  /*   << "y1+" << '\t' << P1_plus_reco[1] << '\n' */
  /*   << "x2+" << '\t' << P2_plus_reco[0] << '\n' */
  /*   << "y2+" << '\t' << P2_plus_reco[1] << '\n' */
  /*   << "x1-" << '\t' << P1_min_reco[0] << '\n' */
  /*   << "y1-" << '\t' << P1_min_reco[1] << '\n' */
  /*   << "x2-" << '\t' << P2_min_reco[0] << '\n' */
  /*   << "y2-" << '\t' << P2_min_reco[1] << '\n' */
  /*   << "z  " << '\t' << z_reco << '\n' */
  /*   << "z_CDA" << '\t' << CDA_mean() << '\n'; */

  delete min;


}

 TMatrixD GetConstraintDerivativeMatrix(Double_t *measures_guess) {
   Double_t x1_p = measures_guess[0];
   Double_t y1_p = measures_guess[1];
   Double_t x1_m = measures_guess[2];
   Double_t y1_m = measures_guess[3];
   Double_t x2_p = measures_guess[4];
   Double_t y2_p = measures_guess[5];
   Double_t x2_m = measures_guess[6];
   Double_t y2_m = measures_guess[7];
   
   TMatrixD D(8,3);

   //Filling with zeroes
   for (int i = 0; i < 8; i++) {
     for (int j = 0; j < 3; j++) {
       D(i,j) = 0;
     }
   }

   //Filling non-zero elements
   D(0,0) = -y2_p;
   D(1,0) = x2_p;
   D(4,0) = y1_p;
   D(5,0) = -x1_p;

   D(2,1) = -y2_m;
   D(3,1) = x2_m;
   D(6,1) = y1_m;
   D(7,1) = -x1_m;

   D(0,2) = x2_m - x1_m;
   D(2,2) = x1_p - x2_p;
   D(4,2) = -x1_m;
   D(6,2) = x1_p;

   /* D(0,0) = 0; */
   /* D(1,0) = 1; */
   /* D(4,0) = 0; */
   /* D(5,0) = 0; */

   /* D(2,1) = 0; */
   /* D(3,1) = 1; */
   /* D(6,1) = 0; */
   /* D(7,1) = 0; */

   /* D(0,2) = 1; */
   /* D(2,2) = 0; */
   /* D(4,2) = 0; */
   /* D(6,2) = 0; */

   return D;
 }

 TMatrixD GetVarianceMatrix() {
   TMatrixD G_inv(8,8);

   for (int i = 0; i < 8; i++) {
     for (int j = 0; j < 8; j++) {
       G_inv(i,j) = 0;
     }
   }

   for (int i = 0; i < 8; i++) {
     G_inv(i,i) = (sigma*sigma);
   }
   
   return G_inv;
 }

 // R = (D_t G^-1 D)^-1
 TMatrixDSym GetRMatrix(TMatrixD D) {

   TMatrixD D_t(3,8);
   D_t.Transpose(D);
   
   TMatrixDSym R(3);

   R.TMult(D);

  
   /* cout << "R: " << endl;  */
   /* cout << "    Det: " << R.Determinant() << endl; */
   
   /* R.Print(); */
 R.InvertFast();
   return R;
 }

 TMatrixD x_Guess(Double_t* measures_guess) {
   
   TMatrixD x_G(8,1);
   
   for (int i = 0; i < 8; i++) {
     x_G(i,0) = measures_guess[i]-measures[i];
   }

   return x_G;
 }
   


 TMatrixD x_Vector (TMatrixD D, TMatrixD x_G, TMatrixD Constr_Vector) {
   
   TMatrixD x_V(8,1);
   TMatrixD R(3,3);
   TMatrixD D_t(3,8);
   D_t.Transpose(D);

   R = GetRMatrix(D);

   TMatrixD A_1(8,3);
   A_1 = D;
   A_1 *= R;
   
   TMatrixD A_2(3,1); //D_t x_G
   A_2.Mult(D_t,x_G);
   
   A_2 -= Constr_Vector;

   x_V.Mult(A_1,A_2);
   
   return x_V;
 }  
   
 void Iterate(Double_t old_coord[], Double_t new_coord[]) {
  
   TMatrixD x_G(8,1);
   TMatrixD x_V(8,1);
   TMatrixD D(8,3);
   /* TMatrixD G_inv(8,8); */
   TMatrixD constr_V(3,1);


   if (isFirstIteration) {
     for (int i = 0; i < 8; i++) {
       old_coord[i] = measures[i];}
     isFirstIteration = false;
   }

   D = GetConstraintDerivativeMatrix(old_coord);
   /* cout << "D: " << endl; */
   /* D.Print(); */
   /* cout << "G_inv: " << endl; */
   /* G_inv = GetVarianceMatrix(); */
   /* G_inv.Print() */;
   
   x_G = x_Guess(old_coord);
   /* cout << "XG" << endl; */
   /* x_G.Print(); */
   /* cout << "FINE XG" << endl; */
   
   constr_V = GetConstraintsVector(old_coord);
   /* cout << "Constr_V: " << endl; */
   /* constr_V.Print(); */
   x_V = x_Vector(D, x_G, constr_V);
   /* cout << "x_V: " << endl; */
   /* x_V.Print(); */
   for (int i = 0; i < 8; i++) {
     new_coord[i] = old_coord[i] + 0.5*(x_V(i,0) - x_G(i,0));
   }


 }
   
 Double_t par_constr[5];

 void IterationCycle(Double_t final_coord[]) {
   
   Double_t old_pars[5] = {0,0,0,0,0}, new_pars[5] = {0,0,0,0,0}, final_pars[5] = {0,0,0,0,0};
   Double_t old_z, new_z;
   Double_t old_coord[8];
   Double_t new_coord[8];
   iteration_no = 0;
   /* ofstream fout("chi2_iter.dat"); */
   /* TFile* file_out = new TFile("iter.root", "recreate"); */
   /* TH1F* chi2_iter_hist = new TH1F("chi2_iter_hist", "Chi2 Iter Hist", 1000,0,1000); */
     /* chi2_iter_hist->SetBinContent(1,chi2_all(measures)); */

 

     for (int j = 0; j < 5; j++) par_constr[j] = 0;
     par_constr[0] = measures[0];
     par_constr[2] = measures[4];
     /* fout << chi2_all(measures) << '\t'; */
     /* fout << z_k(par_constr) << '\t'; */
     /* fout << measures[0] << '\t' << measures[1] << '\t' << measures[3] << endl; */
     /* chi2_iter_hist->SetBinContent(1,chi2_all(measures)); */
     

     do {
     Iterate(old_coord, new_coord);
     iteration_no++;

     old_pars[0] = old_coord[0];
     old_pars[2] = old_coord[4];

     new_pars[0] = new_coord[0];
     new_pars[2] = new_coord[4];

     old_z = z_k(old_pars);
     new_z = z_k(new_pars);
     
     //    cout << "It. no: " << iteration_no << "\t" << new_z << endl;

     for (int j = 0; j < 5; j++) par_constr[j] = 0;
     par_constr[0] = new_coord[0];
     par_constr[2] = new_coord[4];
     /* fout << chi2_all(new_coord) << '\t'; */
     /* fout << z_k(par_constr) << '\t'; */
     /* fout << new_coord[0] << '\t' << new_coord[1] << '\t' << new_coord[3] << endl; */
     /* chi2_iter_hist->SetBinContent(i+2,chi2_all(new_coord)); */

     
     for (int j = 0; j < 8; j++) old_coord[j] = new_coord[j];
     }
     
     while (TMath::Abs(new_z - old_z) > 1E-5 && iteration_no < 5000);
     
     //   cout << "-------------------------------------------" << endl;

     if (iteration_no == 5000) {cout << "MAX iteration number exceeded: decay number " << dec_no << endl; isValid = false;
      }
     
   
   for (int i = 0; i < 8; i++) {
     final_coord[i] = new_coord[i];
   }

   final_pars[0] = final_coord[0];
   final_pars[2] = final_coord[4];
   
   if (iteration_no == 5000) cout << "z: " << z_k(final_pars) << endl;


   /* chi2_iter_hist->Write(); */
   /* file_out->Close(); */
   /* fout.close(); */
 }

};

int event::dec_no = 0;
int event::ev_no = 0;
Double_t event::K_z = 0;
Double_t event::K_p = 0; 
Double_t event::pi_plus_modp = 0; 
Double_t event::pi_plus_theta = 0; 
Double_t event::pi_plus_phi = 0;
Double_t event::pi_min_modp = 0;
Double_t event::pi_min_theta = 0; 
Double_t event::pi_min_phi = 0;
Double_t event::x1_plus = 0; 
Double_t event::y1_plus = 0; 
Double_t event::z1_plus = 0; 
Double_t event::x1_min = 0; 
Double_t event::y1_min = 0; 
Double_t event::z1_min = 0;
Double_t event::x2_plus = 0; 
Double_t event::y2_plus = 0; 
Double_t event::z2_plus = 0; 
Double_t event::x2_min = 0; 
Double_t event::y2_min = 0; 
Double_t event::z2_min = 0;
Double_t event::x1_plus_k = 0; 
Double_t event::y1_plus_k = 0; 
Double_t event::z1_plus_k = 0; 
Double_t event::x1_min_k = 0; 
Double_t event::y1_min_k = 0; 
Double_t event::z1_min_k = 0;
Double_t event::x2_plus_k = 0; 
Double_t event::y2_plus_k = 0; 
Double_t event::z2_plus_k = 0; 
Double_t event::x2_min_k = 0; 
Double_t event::y2_min_k = 0; 
Double_t event::z2_min_k = 0;
Double_t event::K_z_k = 0;

Double_t event::sigma = 0.001;

#endif 

    
