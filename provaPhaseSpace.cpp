#include "FGenPhaseSpace.h"
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TH2F.h>

void PhaseSpace() {

   TLorentzVector target(0.0, 0.0, 0.0, 0.938);
   TLorentzVector beam(0.0, 0.0, .65, .65);
   TLorentzVector W = beam + target;

   //(Momentum, Energy units are Gev/C, GeV)
   Double_t masses[3] = { 0.938, 0.139, 0.139} ;
   /*   Double_t weight;
   TGenPhaseSpace event;
   event.SetDecay(W, 3, masses); */

   FGenPhaseSpace event;
   event.Initialize(W, 3, masses);


   TH2F *h2 = new TH2F("h2","h2", 50,1.1,1.8, 50,1.1,1.8);

   for (Int_t n=0;n<100000;n++) {
     event.GenerateEvent();
     
     //     weight = event.Generate();
     TLorentzVector *pProton = event.GetDecay(0);
     
     TLorentzVector *pPip    = event.GetDecay(1);
     TLorentzVector *pPim    = event.GetDecay(2);
     
     TLorentzVector pPPip = *pProton + *pPip;
     TLorentzVector pPPim = *pProton + *pPim;
     
     h2->Fill(pPPip.M2() ,pPPim.M2());
   }
  
   h2->Draw();
   
}
