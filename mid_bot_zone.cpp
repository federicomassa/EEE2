//LA FINE DI OGNI RIGA DEVE CONTENERE UNO SPAZIO BIANCO
#include <TF1.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <iostream>
#include "fit3.cpp"  //contiene le classi point e triplet
#include <fstream>
#include "cell.h"
#include <string>
#include <sstream>

using namespace std;

int GetHitCount(string str){ //Calcola il numero di hits per evento
  int s = 0;
  for (unsigned int i = 0; i < str.size()+1; i++) {
    if(!isalnum(str[i]) && str[i] != '.' && str[i] != '-') //se trova un carattere non alfanumerico o . allora è la fine di un blocco
      s = s+1;
  }
  return (s-9)/3;
}; // 9 blocchi prima dell'inizio delle 3 coordinate

void mid_bot_tracks(){
  string test;
  double Lx = 82, Ly = 158;
  double x1, y1, z1;
  const int umax = 24, vmax = 2;
  int count2[umax][vmax];
  int count3[umax][vmax];
  double eff[umax][vmax];
  int cellx = 0, celly = 0;
  int check[umax][vmax];
  TH2F* eff_zone = new TH2F("eff_zone","Efficienza a zone, camera 1",96,-39.2917,39.2917,2,-79,79);
  for (int u = 0; u < umax; u++) {
    for (int v = 0; v < vmax; v++) {
      count2[u][v] = 0;
      count3[u][v] = 0;
      eff[u][v] = 0;
      check[u][v] = 0;
    }}

  int evnum = 0;
  // double *besty = new double[3];
  // double parameter = 0;
   ifstream run("../EEEData/CORR_EEE_Prova_SCINT9500__20140603_205538.txt"); //INPUT FILE  
  TFile rfile("top_zone9500.root","RECREATE");
  triplet n1;
  int j = 0;
  double number = 0;
  string line;
  string entry = "";
  int ch1 = 0;
  int ch2 = 0;
  int ch3 = 0;
  int linecount = 0;
  int k = -1;
  do{
    getline(run,line);
  }
  while (line.substr(11,5) != "EVENT");
  // for (int n = 0; n <= 108;n++) {getline(run,line);}

  //   for (k = 0; k <= 10000; k++){

   do  {k+= 1;// cout << "INIZIO DEL DO" << endl;
	 //  if (m%5000 == 0) cout << m << " eventi analizzati..." << endl;
     if (line.find("EVENT") > 15 ) {getline(run,line);continue;}//Durante il run compaiono righe non di evento, se non lo trova find restituisce un numero molto grande
    ch1 = 0;
    ch2 = 0;
    ch3 = 0;

    linecount = GetHitCount(line);

    if (linecount == 0){getline(run,line); continue;}

       point* hit = new point[linecount];//alloca la memoria per tutti i punti dell'evento
       //   cout << point::n << endl;
  // cout << GetHitCount(line) << endl;
  for (unsigned int i = 0; i < line.size()+1; i++) {
    if(!isalnum(line[i]) && line[i] != '.' && line[i] != '-') {
      j = j+1;
      if (entry != "" && j == 4) {stringstream(entry) >> evnum; cout << "EVENTO NUMERO: " << evnum << endl;} 
      if (entry != "" && j >= 9) {
	stringstream(entry) >> number;
	//    	cout << number << endl;
	hit[int(floor((double(j)-9)/3))].SetValue(j%3,number);
	//	    cout << "OK" << endl;
	if (j%3 == 2) {
	  switch (hit[int(floor((double(j)-9)/3))].GetChNumber()) {
	  case(1):
	    ch1 += 1;
	    break;
	  case(2):
	    ch2 += 1;
	    break;
	  case(3):
	    ch3 += 1;
	    break;
	  };}
      } //divide l'evento nei vari punti
      entry = "";}
      else entry = entry + line.substr(i,1);
  } //fine ciclo linea

  

	  
  //   cout << "DOPO FOR" << endl;
 //Controllo la bontà dell'evento: ha almeno un hit per camera (top_mid)?
  if (ch3 == 0 || ch2 == 0){/*cout << "EVENTO NON BUONO" << endl;*/ 
    j = 0; 
    delete [] hit;
    // cout << point::n << endl;
    getline(run,line);
    entry = ""; 

    continue;} //Evento non buono: prossimo evento
  
  else {//evento con almeno un hit per camera (top_mid)
  

    if (ch1==0) { //funziona solo con file corretti
      for (int a = 0; a < ch2; a++) {
	for (int b = 0; b < ch3; b++) {
	  z1 = zch1;
	  x1 = ZXGetIntercept(hit[a],hit[b+ch2]) + ZXGetSlope(hit[a],hit[b+ch2])*z1;
	  y1 = (z1-YZGetIntercept(hit[a],hit[b+ch2]))/YZGetSlope(hit[a],hit[b+ch2]);
	

	  if ((sqr(x1) < sqr(82.0) && (sqr(y1) < sqr(158.0)))) {
	      cellx = xcell(x1,Lx,umax);
	      celly = ycell(y1,Ly,vmax);
	      if (check[cellx][celly] == 0) count2[cellx][celly] += 1;
	      check[cellx][celly] = 1;
	  }
	  
	      
	}
      }
      
    }
  


    if(ch1 >= 1){
       //    for (int kk = 0; kk < 3; kk++) {
    //	cout << hit[kk].x << endl;
    //	cout << hit[kk].y << endl;
    //	cout << hit[kk].z << endl;}
   for (int a = 0; a < ch1; a++) {
     for (int b = ch1; b < ch1+ch2; b++) {
       for (int c = ch1+ch2; c < ch1+ch2+ch3; c++) {
   	 	 n1.SetPoints(hit[a],hit[b],hit[c]); //considero tutte le combinazioni di triplette
		 cellx = xcell(hit[a].x,Lx,umax);
		 celly = ycell(hit[a].y,Ly,vmax);
		 if (check[cellx][celly] == 0){ count3[cellx][celly] += 1;}
		 check[cellx][celly] = 1;
   	
       }}}
    }
  

    for (int u = 0; u < umax; u++) {
      for (int v = 0; v < vmax; v++) {
	check[u][v] = 0;
      }
    }

  j = 0;
  delete[] hit;
  getline(run,line);
  entry = "";
  continue;
  }

    //vedere se la traccia entra nella terza camera o meno!
  j = 0;

  delete[] hit;
  // cout << point::n << endl;
  getline(run,line);
  entry = "";
   
   }
 // cout << "FINE DEL DO" << endl;



     	  	  while (!run.eof());

   for(int u = 0; u < umax; u++) {
     for(int v = 0; v < vmax; v++) {
       eff[u][v] = double(count3[u][v])/double(count3[u][v]+count2[u][v]);
       eff_zone->SetBinContent((u+1)*4,v+1,eff[u][v]);
     }
   }


   for (int u = 0; u < umax; u++) {
     for (int v = 0; v < vmax; v++) {
       if (v == vmax -1) cout << eff[u][v] << endl;
       else cout <<  eff[u][v] << '\t';
     }
   }
   eff_zone->SetDrawOption("COLZ");
   eff_zone->GetZaxis()->SetRangeUser(0.65,1);
   eff_zone->Write();
   rfile.Close();
     run.close();
    
   
}
  
