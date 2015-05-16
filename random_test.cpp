//TEST DI CASUALITA' DEL GENERATORE TRANDOM3

#include <TMath.h> //libreria ROOT matematica
#include <TRandom3.h> //libreria ROOT generatore pseudocasuale
#include <TH1F.h> //libreria ROOT istogrammi 1D
#include <TFile.h> //libreria ROOT per scrittura file
#include <TCanvas.h> //libreria ROOT per creazione di canvas ("tele" su cui vengono rappresentati i grafici)
#include <iostream>


//invece che usare una classe avrei potuto dichiarare una variabile globale rndgenerator
//ma questa struttura è utile per il generatore di eventi dello spettrometro
class ev_generator {

private: //i membri private sono visibili solo dall'interno della classe o dalle classi "amiche". Evita che l'utente possa accedere direttamente al generatore
  TRandom3 rndgenerator; //random number generator

public: //tutti i membri dichiarati come public sono visibili dall'esterno della classe



//funzione che restituisce la coppia di punti x,y tra 0 e 1
//i & servono per passare gli argomenti come "referenza" e non come valore, in modo tale che siano modificabili dall'interno della funzione
  void generate(double &x, double &y) { 

    x = rndgenerator.Uniform();
    y = rndgenerator.Uniform();
 
  }
    

};


void random_test() {
  
  ev_generator event; //oggetto di tipo ev_generator 
  int acc = 0; //conta il numero di accettati
  const int N = 2000; //numero di esperimenti generati
  const int imax = 100000; //numero di coppie generate per ogni esperimento
  double x, y; //valori x,y generati
  double pi[N]; //array che contiene tutti gli N valori di pi simulati
  double mean_pi = 0; //media di pi negli N esperimenti
  double var_pi = 0, stddev_pi; //var e std_dev di pi negli N esperimenti

  TFile* f_out = new TFile("random_test.root","RECREATE"); //file root in uscita, da visualizzare su TBrowser. RECREATE sovrascrive il file se già esistente
  TH1F* pi_hist = new TH1F("pi_hist","Istogramma valori di pi simulati;pi;#",20, 3.11,3.18); //dichiarazione istogramma: nome, titolo;asse_x;asse_y, numero_bin, estremo sx, estremo dx

  //genero imax eventi e calcolo pi quando i punti sono sotto il cerchio di raggio 1
  //calcolo N volte pi e 
  for (int n = 0; n < N; n++) {
    for (int i = 0; i < imax; i++) {
      event.generate(x,y); //generazione evento
      if (x*x+y*y < 1) acc += 1;  //accetto solo se sotto il cerchio
    }
    pi[n] = 4*double(acc)/double(imax); //  PI/4 = eventi_accettati/eventi_totali
    acc = 0; //reset
    
    pi_hist->Fill(pi[n]);
  }


  //Calcolo la media
  for (int i = 0; i < N; i++) {
    mean_pi += pi[i]/N;
  }

  //Calcolo la varianza
  for (int i = 0; i < N; i++) {
    
    var_pi += TMath::Power((pi[i]-mean_pi),2)/(N-1); //funzione potenza

  }

  var_pi = var_pi/N; //var sulla media

  stddev_pi = TMath::Sqrt(var_pi); //funzione radice

  std::cout << "PI: " << mean_pi << " STD_DEV: " << stddev_pi << std::endl; //stampa il risultato su terminale
  
  TCanvas* canv = new TCanvas(); //nuovo canvas senza titolo
  canv->cd(); //imposta come current canvas, non necessario ma toglie warning di compilazione
  pi_hist->DrawCopy(); //disegna su canvas canv
  pi_hist->Write(); //scrive su f_out
  
  f_out->Close(); //chiude file
}
