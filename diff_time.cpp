//considero sorgente con rate 1 Hz e guardo i primi 10 s

#include <TRandom3.h>
#include <TH1F.h>
#include <algorithm>
#include <TFile.h>

using namespace std;

void diff_time() {

  double r[50000];
const int imax = 50000;
 TRandom3 rndgen;
 TH1F* dist = new TH1F("dist","Distribuzione degli intervalli di tempo tra due eventi consecutivi; Intervallo di tempo (s); #",100,0,10);
 TFile* out = new TFile("dist_random.root","RECREATE");
 for (int i = 0; i < imax; i++) {
r[i] = rndgen.Uniform(50000);
 }

sort(r,r+imax);

 for (int i = 0; i < imax-1; i++) {
   dist->Fill(r[i+1]-r[i]);
 }

 dist->Write();
 out->Close();
}


 

		       
