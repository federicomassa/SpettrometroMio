#include <TMath.h>
#include <TMinuit.h>
#include <iostream>

using namespace std;

double P1_plus[3], P2_plus[3], P1_min[3], P2_min[3];

double y2_plus(Double_t *par) {
  return par[2]/par[0]*par[1];
}

double x2_min(Double_t *par) {
  return par[3]/par[0]*par[2];
}

double y2_min(Double_t *par) {
  return par[2]/par[0]*par[4];
}

double z_dec(Double_t *par) {
  return P1_plus[2] - (P2_plus[2]-P1_plus[2])/(par[2]-par[0])*par[0];
}


void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  
  ///////// Parameters: x1+, y1+, x2+, x1-, y1-. The rest has been removed due to constraints.
  double chi;
  double sigma = 0.001;
  

  chi = 
    TMath::Power((par[0] - P1_plus[0])/sigma,2) + 
    TMath::Power((par[1] - P1_plus[1])/sigma,2) + 
    TMath::Power((par[2] - P2_plus[0])/sigma,2) + 
    TMath::Power((y2_plus(par) - P2_plus[1])/sigma,2) +
    TMath::Power((par[3] - P1_min[0])/sigma,2) + 
    TMath::Power((par[4] - P1_min[1])/sigma,2) + 
    TMath::Power((x2_min(par) - P2_min[0])/sigma,2) + 
    TMath::Power((y2_min(par) - P2_min[1])/sigma,2); 
  
  f = chi;

}


void kinematic_z_reco() {

  double P1_plus_reco[3], P2_plus_reco[3], P1_min_reco[3], P2_min_reco[3];
  double final_pars[5], final_pars_err[5];
  double z_reco;

  P1_plus[0] = 0.08372;	
  P1_plus[1] = -0.02893;	
  P1_plus[2] = 25.00000;	
  P2_plus[0] = 0.13810;	
  P2_plus[1] = -0.04734;	
  P2_plus[2] = 35.00000;	

  P1_min[0] = -0.03807;
  P1_min[1] = 0.01335;
  P1_min[2] = 25.00000;	
  P2_min[0] = -0.06537;	
  P2_min[1] = 0.02213;
  P2_min[2] = 35.00000;

  TMinuit* gMinuit = new TMinuit(5);
  gMinuit->SetFCN(chi2);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  Double_t vstart[5];
  vstart[0] = P1_plus[0];
  vstart[1] = P1_plus[1];
  vstart[2] = P2_plus[0];
  vstart[3] = P1_min[0];
  vstart[4] = P1_min[1];
  
  Double_t step[5] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
  
  gMinuit->mnparm(0, "x1+", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "y1+", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "x2+", vstart[2], step[2], 0,0,ierflg);
  gMinuit->mnparm(3, "x1-", vstart[3], step[3], 0,0,ierflg);
  gMinuit->mnparm(4, "y1-", vstart[4], step[4], 0,0,ierflg);

  arglist[0] = 500;
  arglist[1] = 1;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

 
  
  gMinuit->GetParameter(0, final_pars[0], final_pars_err[0]); 
  gMinuit->GetParameter(1, final_pars[1], final_pars_err[1]);
  gMinuit->GetParameter(2, final_pars[2], final_pars_err[2]);
  gMinuit->GetParameter(3, final_pars[3], final_pars_err[3]);
  gMinuit->GetParameter(4, final_pars[4], final_pars_err[4]);

  P1_plus_reco[0] = final_pars[0];
  P1_plus_reco[1] = final_pars[1];
  P2_plus_reco[0] = final_pars[2];
  P2_plus_reco[1] = y2_plus(final_pars);
  
  P1_min_reco[0] = final_pars[3];
  P1_min_reco[1] = final_pars[4];
  P2_min_reco[0] = x2_min(final_pars);
  P2_min_reco[1] = y2_min(final_pars);

  z_reco = z_dec(final_pars);

  cout << "FINAL RESULTS " << '\n' 
    << "x1+" << '\t' << P1_plus_reco[0] << '\n'
    << "y1+" << '\t' << P1_plus_reco[1] << '\n'
    << "x2+" << '\t' << P2_plus_reco[0] << '\n'
    << "y2+" << '\t' << P2_plus_reco[1] << '\n'
    << "x1-" << '\t' << P1_min_reco[0] << '\n'
    << "y1-" << '\t' << P1_min_reco[1] << '\n'
    << "x2-" << '\t' << P2_min_reco[0] << '\n'
    << "y2-" << '\t' << P2_min_reco[1] << '\n'
    << "z  " << '\t' << z_reco << '\n'
    << "z_CDA" << '\t' << CDA_mean() << '\n';
  

}
