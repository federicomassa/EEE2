#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH2F.h>
#include "fit3.cpp"
#include "cell.h"

//lati lunghi = 48.5 cm, lati corti = 40 cm, spessore = 1.2 cm, distanze relative 10,2 cm (scintillatori per flusso raggi cosmici)
double xmax = 82; //cm //confermate da EEEpaper
double ymax = 158; //cm
double D12 = 53.2;// cm //misure prese da mail 
double s_D12 = 0.00001;
double D23 = 52.8; //cm
double s_D23 = 0.00001;

int imax = 1000000;
double strd = 82.0/24.0; //unità arbitrarie, distanza tra due strip considerate filiformi ---- CORRETTA -----
double stheta;
double sphi;
using namespace std;


void stripaccdiscreta(){

   double eps1[24][2] = {{0.730964,0.680162},{0.877301,0.847926},{0.908497,0.875969},{0.883721,0.851064},{0.935484,0.854545},{0.866359,0.871479},{0.931973,0.85559},{0.873646,0.813137},{0,0},{0.880524,0.823642},{0.9046,0.872165},{0.879339,0.824338},{0.912249,0.875823},{0.896277,0.820717},{0.92673,0.851738},{0.894488,0.830769},{0.903361,0.883469},{0.813609,0.830153},{0.837838,0.866412}/*questa a sinistra in teoria è rotta*/,{0.837838,0.866412},{0.887324,0.921397},{0.843137,0.907489},{0,0},{0.705584,0.70852}};
   double eps2[24][2] = {{0.979167,1},{0.984375,0.982759},{0.973913,0.969925},{0.946809,0.967914},{0.958678,0.934307},{0.912162,0.895735},{0.854749,0.865823},{0.794595,0.890566},{0.741514,0.748727},{0,0},{0.777778,0.757202},{0.778027,0.822835},{0.760606,0.850746},{0.85814,0.830657},{0.865714,0.84854},{0.795302,0.863636},{0.797414,0.913043},{0.819767,0.9},{0.918699,0.944785},{0.943925,0.993464},{0.991453,0.978102},{0.989899,0.992},{0.981132,1},{0.98913,1}};
   double eps3[24][2] = {{0.756757,0.78125},{0.962025,0.878505},{0.976471,0.936937},{0.938144,0.903509},{0.952381,0.959184},{0.955307,0.920732},{0.962733,0.911765},{0.959288,0.91369},{0.935691,0.896203},{0.917763,0.886889},{0.938838,0.89589},{0.902174,0.898667},{0.941781,0.918782},{0.909639,0.894253},{0.890566,0.892857},{0.835821,0.793939},{0,0},{0.790698,0.771144},{0,0},{0.797101,0.850932},{0,0},{0.905882,0.872881},{0.941176,0.900901},{0.705882,0.769231}}; 

  // double eps1[24][2];
  // double eps2[24][2];
  // double eps3[24][2];

  // for (int n1 = 0; n1 < 24; n1++) {
  //   for (int n2 = 0; n2 < 2; n2++) {
  //     eps1[n1][n2] = 1;
  //     eps2[n1][n2] = 1;
  //     eps3[n1][n2] = 1;
  //   }
  // }


 int ny1 = 0, ny2 = 0, ny3 = 0;
 //  double rD12, rD23;
  point p1,p2,p3;
triplet n1;
	TH1F* htheta = new TH1F("dis_acctheta","Distribuzione Theta accettati", 50, 0, 90);
        TH1F* hrestheta = new TH1F("dis_restheta","Distribuzione Theta-Theta_quantizzato accettati", 50, -45, -45);
	TH1F* hstheta = new TH1F("discr_dis_acctheta","Distribuzione Theta accettati (strip)", 50, 0, 90);

	TH2F* dthvsth = new TH2F("dthvsth","Stheta-Theta vs. Theta Correlation;stheta-theta(°);theta(°)",400,-4,4,60,0,60);
	TH2F* dphvsph = new TH2F("dphvsph","Sphi-phi vs. phi Correlation;sphi-phi(°);phi(°)",500,-100,100,360,0,360);
	TH2F* dphcritvsth = new TH2F("dphcritvsth","Sphi-phi vs. Theta Correlation, sPhi=90,270;theta(°);sphi-phi(°)",500,0,40,500,-100,100);

        TH1F* hresphi = new TH1F("dis_resphi","Distribuzione Phi-Phi_quantizzato accettati", 40, -20, 20);
	TH1F* hphi = new TH1F("dis_accphi", "Distribuzione Phi accettati", 90, 0, 360);	
	TH1F* hsphi = new TH1F("discr_dis_accphi", "Distribuzione Phi accettati (strip)", 90, 0, 360);
	TH2F* phitheta = new TH2F("phitheta","Phi-Theta Correlation;phi(°);theta(°)",90,0,360,50,0,90);
	TH2F* sphistheta = new TH2F("sphistheta","sPhi-sTheta Correlation;sphi(°);stheta(°)",90,0,360,50,0,90);
	TFile rfile("accettanza_missing_strips_zone.root","RECREATE");
	TRandom3 rndgen;
  ofstream* dtheta = new ofstream("theta.dat");
  ofstream* acctheta = new ofstream("acctheta.dat");
  ofstream* accphi = new ofstream("accphi.dat");
  //efficienze
  //  float eps1 [3][3] = {{0.612,0.662,0.762},{0.425,0.470,0.612},{0.574,0.618,0.722}};
  double theta, W1,W2,W3, phi,x1,x2,x3,y1,y2,y3,z1,z2 = 0,z3,xs1,xs2,xs3;
  int j = 0;
  int ns1, ns2,ns3;
 
  for (int i = 1; i <= imax;) {
    theta = rndgen.Uniform(2*atan(1.)); 
    W1 = rndgen.Uniform(1);
    if (sin(theta)*pow(cos(theta),2) > W1) {
        i+= 1;
      *dtheta << theta << endl;
      
      if (i%5000 == 0) cout << double(i)/double(imax)*100 << "%" << endl;
      phi = rndgen.Uniform(8*atan(1.));
      x2 = (rndgen.Uniform(xmax)-xmax/2.);
      y2 = (rndgen.Uniform(ymax)-ymax/2.);
      x1 = tan(theta)*cos(phi)*D12 + x2;
      y1 = tan(theta)*sin(phi)*D12 + y2;
      z1 = D12;
      x3 = -tan(theta)*cos(phi)*D23+x2;
      y3 = -tan(theta)*sin(phi)*D23 + y2;
      z3 = -D23;
      W1 = rndgen.Uniform(1);
      W2 = rndgen.Uniform(1);
      W3 = rndgen.Uniform(1);

      ns1=xcell(x1,xmax,24);
      ns2=xcell(x2,xmax,24);
      ns3=xcell(x3,xmax,24);

      ny1 = ycell(y1,ymax,2);
      ny2 = ycell(y2,ymax,2);
      ny3 = ycell(y3,ymax,2);


	xs1=xstrip(x1);
	xs2=xstrip(x2);
	xs3=xstrip(x3);
      
	if ((pow(x1,2) <= pow(xmax,2)/4) && (pow(y1,2) <= pow(ymax,2)/4) && (pow(x3,2) <= pow(xmax,2)/4) && (pow(y3,2) <= pow(ymax,2)/4)) {
	
    
	     if(W1 <= eps1[ns1][ny1] && W2 <= eps2[ns2][ny2] && W3 <= eps3[ns3][ny3]) {
	    

	


	p1.SetValues(xs1,y1,z1);
	p2.SetValues(xs2,y2,z2);
	p3.SetValues(xs3,y3,z3);

	
	n1.SetPoints(p1,p2,p3);
	
	n1.Fit();

	stheta=n1.GetTheta();

	sphi=n1.GetPhi();
	  
	if (ns1 !=-3 && ns1 !=10 && ns2 !=-2 && ns3 != 8 && ns3 != 6 && ns3 !=4) {
	
	  *acctheta << theta << endl;
	  *accphi << phi << endl;

	  htheta->Fill(theta*180/3.14159);
	  hstheta->Fill(stheta*180/3.14159);
	  hrestheta->Fill((stheta-theta)*180/3.1415);
	  dthvsth->Fill((stheta-theta)*180/3.1415,theta*180/3.14159);

	  hphi->Fill(phi*180/3.14159);
	  hsphi->Fill(sphi*180/3.1415);
	  hresphi->Fill((sphi-phi)*180/3.14159);
	  dphvsph->Fill((sphi-phi)*180/3.14159,phi*180/3.14159);
	  phitheta->Fill(phi*180/3.14159,theta*180/3.14159);
	  sphistheta->Fill(sphi*180/3.14159,stheta*180/3.14159);
	  j+=1;
	}
		}
		 }
  }

   }		


  htheta->Scale(21976/double(j));
  hphi->Scale(21976/double(j));
  hstheta->Scale(21976/double(j));
  hsphi->Scale(21976/double(j));

  phitheta->Write();
  sphistheta->Write();

  htheta->Write();
  hstheta->Write();
  hrestheta->Write();

  hphi->Write();
  hsphi->Write();
  hresphi->Write();

  dthvsth->Write();
  dphvsph->Write();

  dphcritvsth->Write();

  rfile.Close();
  dtheta->close();
  acctheta->close();
  accphi->close();
  cout << "Generati: " << imax << endl;
  cout << "Accettati: " << j << endl;
  cout << "Accettanza: " << double(j)/double(imax) << endl;
}
