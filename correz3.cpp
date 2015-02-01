//corregge il file prendendo come estremi della strip quelli che danno due zeri di seguito
//LA FINE DI OGNI RIGA DEVE CONTENERE UNO SPAZIO BIANCO
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <iostream>
#include "fit.cpp"  //contiene le classi point e triplet
#include <fstream>

#include <string>
#include <sstream>

using namespace std;

int GetHitCount(string str){ //Calcola il numero di hits per evento
  int s = 0;
  for (unsigned int i = 0; i < str.size(); i++) {
    if(!isalnum(str[i]) && str[i] != '.' && str[i] != '-') //se trova un carattere non alfanumerico o . allora è la fine di un blocco
      s = s+1;
  }
  return (s-9)/3;
}; // 9 blocchi prima dell'inizio delle 3 coordinate



void correz(){
  double temp = 0;
  int stripnum = -1;
  double corrx = 0;
  double corry1[24], corry2[24], corry3[24];

  for (int i = 0; i < 24; i++) {corry1[i] = 0; corry2[i] = 0; corry3[i] = 0;}

  int chnum = 0;
  char dest1[80] = "../EEEData/";
  char dest2[80] = "../EEEData/";
  char ofile[80];
  strcpy(ofile,dest2);
   char infile1[80] = "EEE_PISA01TestRun4Telescopes_20140507_014455.txt";
 char infile2[80] = "EEE_PISA01TestRun4Telescopes_20140507_014455.txt";
   // char infile2[80] = "EEE_Prova_topbottom8900__20140530_174533.txt";
   //  char infile2[80] = "EEE_Prova_topmid9000__20140530_181026.txt";
   //   char infile2[80] = "EEE_Prova_midbottom9000__20140530_184644.txt"; //questo è quello che andrà corretto
  strcat(dest1,infile1);
  strcat(dest2,infile2);
  ifstream run(dest1); //INPUT FILE 
  ifstream tobecorr(dest2);
  char crr[80] = "CORR_";
  strcat(ofile,crr);
  strcat(ofile,infile2);
  ofstream crrfile(ofile);
  int check = 0;
  int bc1, bc2, bc3;
  double up1[24],up2[24],up3[24],low1[24],low2[24],low3[24];
  int evnum = 0;
 
  TFile rfile("correzione_andall.root","RECREATE");
  TH2F* disxy1 = new TH2F("disxy1","XY Occupancy: Ch 1", 200,-100,100,800,-800,800);
  TH2F* disxy2 = new TH2F("disxy2","XY Occupancy: Ch 2", 200,-100,100,800,-800,800);
  TH2F* disxy3 = new TH2F("disxy3","XY Occupancy: Ch 3", 200,-100,100,800,-800,800);
  // triplet n1;
  int j = 0;
  double number = 0, numbery = 0, numberx = 0;
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

  //    for (k = 0; k <= 50000; k++){

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
 //Controllo la bontà dell'evento: ha almeno un hit per camera (top_bottom)?
  // if (ch1 < 1 ||ch2 < 1 || ch3 < 1){/*cout << "EVENTO NON BUONO" << endl;*/ 
  //   j = 0; 
  //   delete [] hit;
    
  //   getline(run,line);
  //   entry = ""; 

  //   continue;} //Evento non buono: prossimo evento

  // else {
    //evento con almeno un hit per camera
 
   for (int a = 0; a < (ch1); a++) {
     disxy1->Fill(hit[a].x,hit[a].y);}

   for (int a = ch1; a < ch1+ch2; a++) {
     disxy2->Fill(hit[a].x,hit[a].y);}

   for (int a = ch1+ch2; a < ch1+ch2+ch3; a++) {
     disxy3->Fill(hit[a].x,hit[a].y);}







  j = 0;
  delete[] hit;
  getline(run,line);
  entry = "";
  //Bin per chamber histograms
   

   }



     	  while (!run.eof());
  
 for (int a = 0; a < 24; a++) {
     for (int b = 1; b <= 400; b++) {
       bc1 = disxy1->GetBinContent(41+5*a,b+400);
       if (bc1 == 0) {if(check == 0) temp = double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     up1[a] = temp; temp = 0; check = 0; 
      

   
     for (int b = 1; b <= 400; b++) {
       bc1 = disxy1->GetBinContent(41+5*a,401-b);
       if (bc1 == 0) {if(check == 0) temp = -double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     low1[a] = temp; temp = 0; check = 0;
 }


 for (int a = 0; a < 24; a++) {
     for (int b = 1; b <= 400; b++) {
       bc2 = disxy2->GetBinContent(41+5*a,b+400);
       if (bc2 == 0) {if(check == 0) temp = double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     up2[a] = temp; temp = 0; check = 0; 
      

   
     for (int b = 1; b <= 400; b++) {
       bc2 = disxy2->GetBinContent(41+5*a,401-b);
       if (bc2 == 0) {if(check == 0) temp = -double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     low2[a] = temp; temp = 0; check = 0;
   }



 for (int a = 0; a < 24; a++) {
     for (int b = 1; b <= 400; b++) {
       bc3 = disxy3->GetBinContent(41+5*a,b+400);
       if (bc3 == 0) {if(check == 0) temp = double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     up3[a] = temp; temp = 0; check = 0; 
      

   
     for (int b = 1; b <= 400; b++) {
       bc3 = disxy3->GetBinContent(41+5*a,401-b);
       if (bc3 == 0) {if(check == 0) temp = -double(b)*2; check += 1; if(check == 2) break;}
       else check = 0;
     }
     low3[a] = temp; temp = 0; check = 0;
   }



  //   for (int bb = 0; bb < 24; bb++) {
 // cout << up1[bb]-low1[bb] << endl;
 // cout << "UP: " << up1[bb] << endl;
 // cout << "LOW: " << low1[bb] << endl;
 // cout << '\n';
 // cout << up2[bb]-low2[bb] << endl;
 // cout << "UP: " << up2[bb] << endl;
 // cout << "LOW: " << low2[bb] << endl;
 // cout << '\n';
 // cout << up3[bb]-low3[bb] << endl;
 // cout << "UP: " << up3[bb] << endl;
 // cout << "LOW: " << low3[bb] << endl;
 // cout << '\n' << '\n';
 //  }

  corrx = (82.0-82.0/24.0)/(60.0+55.0);  // 82-82/24 è la nuova distanza tra le strip, 60+55 la vecchia
  for (int i = 0; i < 24; i++) {corry1[i] = up1[i]-low1[i]; corry2[i] = up2[i]-low2[i]; corry3[i] = up3[i] - low3[i];}

  for (int i = 0; i < 24; i++) {
   corry1[i] = 158.0/corry1[i];  //y rescaling
   corry2[i] = 158.0/corry2[i];
   corry3[i] = 158.0/corry3[i];}

   cout << '\n' << "FINITO. ORA RICOMINCIA E CORREGGE" << '\n' << endl;

   ////////////////////////////////////RICOMINCIA DA CAPO E CORREGGI ///////////////////////////////////////
   run.close();
   k = -1;
   j = 0;

   getline(tobecorr,line);

   // for (k = 0; k <= 50000; k++){

   do  {k+= 1;// cout << "INIZIO DEL DO" << endl;
	 //  if (m%5000 == 0) cout << m << " eventi analizzati..." << endl;
     if (line.find("EVENT") > 15 ) {crrfile << line << endl; getline(tobecorr,line);j = 0; continue;}//Durante il run compaiono righe non di evento, se non lo trova find restituisce un numero molto grande
    linecount = GetHitCount(line);
    if (linecount == 0){getline(tobecorr,line); j = 0; continue;}

      point* hit = new point[linecount];//alloca la memoria per tutti i punti dell'evento
       //   cout << point::n << endl;
  // cout << GetHitCount(line) << endl;
  for (unsigned int i = 0; i < line.size()+1; i++) {
    if(!isalnum(line[i]) && line[i] != '.' && line[i] != '-') {
      j += 1;
      if (j < 9) crrfile << entry << ' ';
      if (entry != "" && j == 4) {stringstream(entry) >> evnum; cout << "EVENTO NUMERO: " << evnum << endl;} 
      if (entry != "" && j >= 9) {
	stringstream(entry) >> number;
	//  cout << number << endl;
	hit[int(floor((double(j)-9)/3))].SetValue(j%3,number);
	//	    cout << "OK" << endl;
	if (j%3 == 2) {
    	chnum = hit[int(floor((double(j)-9)/3))].GetChNumber();
	if (chnum == 1) {number = 53.2; crrfile << (numbery - (up1[stripnum]+low1[stripnum])/2)*corry1[stripnum] << ' ';}
	if (chnum == 2) {number = 0; crrfile << (numbery - (up2[stripnum]+low2[stripnum])/2)*corry2[stripnum] << ' ';}
	if (chnum == 3) {number = -52.8; crrfile << (numbery - (up3[stripnum]+low3[stripnum])/2)*corry3[stripnum] << ' ';}
	crrfile << number << ' ';
	} 
	if (j%3 == 0) {numberx = number; stripnum = int((numberx+60.0)/5.0); crrfile << (numberx*corrx+41.0/24.0) << ' ';}
	if (j%3 == 1) numbery = number;
      }
      entry = "";}
      else entry = entry + line.substr(i,1);
  } //fine ciclo linea
     j = 0;
  delete[] hit;
  getline(tobecorr,line);
  entry = "";
  crrfile << endl;}
   while (!tobecorr.eof());
    
   disxy1->Write();
   disxy2->Write();
   disxy3->Write();
     rfile.Close();
     tobecorr.close();
     crrfile.close();
}
  
