//LA FINE DI OGNI RIGA DEVE CONTENERE UNO SPAZIO BIANCO
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <iostream>
#include "fit3.cpp"  //contiene le classi point e triplet con assi scambiati, usando LSM
#include <fstream>

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



void tracks(){  
  bool bestvalid = false;
  bool bestxvert = false;
  string test;
  int besta = 100, bestb = 100, bestc = 100;
  //  double *besty = new double[3];
  // ifstream run("../Data/EEE_Prova_topbottom8900__20140530_174533.txt"); //INPUT FILE
    ifstream run("../EEEData/CORR_EEE_PISA01TestRun4Telescopes_20140507_014455.txt"); //INPUT FILE
  //  int point::n = 0;
    // TFile rfile("Disttopbot.root","RECREATE");
 TFile rfile("Dist_all_trasl.root","RECREATE");
 TH1F* time_dist = new TH1F("time_dist", "Distribuzione dei tempi tra due eventi consecutivi; Tempo (ms);# ", 1000,0,1000);
 TF1* manyhitsline = new TF1("prova","[0]+[1]*x",0,15);
  TH2F* manyhits2 = new TH2F("manyhits","Many-hits events correlation: Ch 2", 30,0,15,100,0,50);
  TH1F* dischi = new TH1F("dischi","Chi2 distribution; chi2; #", 100,0,10);
TH1F* disyxchi = new TH1F("disyxchi","YX Chi2 distribution; chi2; #", 300,0,30);
TH1F* diszxchi = new TH1F("diszxchi","ZX Chi2 distribution; chi2; #", 300,0,30);
TH1F* disyzchi = new TH1F("disyzchi","YZ Chi2 distribution; chi2; #", 300,0,30);
  TH1F* hpc1 = new TH1F("hpc1", "Hit per chamber / Chamber 1; #Hits;# ", 20,0,20);
  TH1F* hpc2 = new TH1F("hpc2", "Hit per chamber / Chamber 2;#Hits;#", 20,0,20);
  TH1F* hpc3 = new TH1F("hpc3", "Hit per chamber /Chamber 3;#Hits;#", 20,0,20);
  TH1F* disprimatheta = new TH1F("disprimatheta","Theta distribution;Theta(deg);#", 50, 0,90);
  TH1F* disprimaphi = new TH1F("disprimaphi","Phi distribution;Phi(deg);#", 90, 0, 360);
  TH2F* disxy1 = new TH2F("disxy1","XY Occupancy: Ch 1", 200,-100,100,100,-400,400);
  TH1F* disx1 = new TH1F("disx1", "X Occupancy: Ch 1", 200,-100,100);
  TH1F* disy1 = new TH1F("disy1", "Y Occupancy: Ch 1", 50,-400,400);
  TH2F* disxy2 = new TH2F("disxy2","XY Occupancy: Ch 2", 200,-100,100,100,-400,400);
  TH1F* disx2 = new TH1F("disx2", "X Occupancy: Ch 2", 200,-100,100);
  TH1F* disy2 = new TH1F("disy2", "Y Occupancy: Ch 2", 50,-400,400);
  TH2F* disxy3 = new TH2F("disxy3","XY Occupancy: Ch 3", 200,-100,100,100,-400,400);
  TH1F* disx3 = new TH1F("disx3", "X Occupancy: Ch 3", 200,-100,100);
  TH1F* disy3 = new TH1F("disy3", "Y Occupancy: Ch 3", 50,-400,400);
  TH1F* disz = new TH1F("disz", "Z Occupancy; Ch Number; #", 15,0,4);
  TH1F* disdist1 = new TH1F("disdist1","Distance distribution, Ch 1;Distance;#",200,0,600);
  TH1F* disdist2 = new TH1F("disdist2","Distance distribution, Ch 2;Distance;#",200,0,600);
  TH1F* disdist3 = new TH1F("disdist3","Distance distribution, Ch 3;Distance;#",200,0,600);
  TH1F* disdistx1 = new TH1F("disdistx1","X Distance distribution, Ch1;X Distance;#",200,-100,100);
  TH1F* disdistx2 = new TH1F("disdistx2","X Distance distribution, Ch2;X Distance;#",200,-100,100);
  TH1F* disdistx3 = new TH1F("disdistx3","X Distance distribution, Ch3;X Distance;#",200,-100,100);
  TH1F* disdisty1 = new TH1F("disdisty1","Y Distance distribution, Ch1;Y Distance;#",800,-400,400);
  TH1F* disdisty2 = new TH1F("disdisty2","Y Distance distribution, Ch2;Y Distance;#",800,-400,400);
  TH1F* disdisty3 = new TH1F("disdisty3","Y Distance distribution, Ch3;Y Distance;#",800,-400,400);
  TH1F* av_dist1 = new TH1F("av_dist1","Average Distance distribution, Ch 1;Distance (cm);#",200,0,600);
  TH1F* av_dist2 = new TH1F("av_dist2","Average Distance distribution, Ch 2;Distance (cm);#",200,0,600);
  TH1F* av_dist3 = new TH1F("av_dist3","Average Distance distribution, Ch 3;Distance (cm);#",200,0,600);
  TH1F* av_distx1 = new TH1F("av_distx1","Average X Distance distribution, Ch 1;Distance (cm);#",200,-100,100);
  TH1F* av_distx2 = new TH1F("av_distx2","Average X Distance distribution, Ch 2;Distance (cm);#",200,-100,100);
  TH1F* av_distx3 = new TH1F("av_distx3","Average X Distance distribution, Ch 3;Distance (cm);#",200,-100,100);
  TH1F* av_disty1 = new TH1F("av_disty1","Average Y Distance distribution, Ch 1;Distance (cm);#",800,-400,400);
  TH1F* av_disty2 = new TH1F("av_disty2","Average Y Distance distribution, Ch 2;Distance (cm);#",800,-400,400);
  TH1F* av_disty3 = new TH1F("av_disty3","Average Y Distance distribution, Ch 3;Distance (cm);#",800,-400,400);
  TH2F* primathetaphi = new TH2F("prima_phi-teta","Phi-Theta Correlation before clustering;Phi(deg);Theta(deg)", 90,0,360,50,0,90);
  TH2F* piccolitheta = new TH2F("piccolitheta","Piccoli theta;Phi(deg);Theta(deg)", 90,0,360,50,0,0.1);
  triplet n1;
  double chi = 1E30, tempchi = 0,yxtempchi = 0, zxtempchi = 0,  yztempchi = 0, theta = 1000, phi = 1000;
  int j = 0;
  double number = 0;
  string line;
  string entry = "";
  int ch1 = 0;
  int ch2 = 0;
  int ch3 = 0;
  int linecount = 0;
  int k = -1;
  int k1 = 0, k2 = 0, k3 = 0;
  double d1 = 0, d2 = 0, d3 = 0, dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0, dx3 = 0, dy3 = 0;
  double seconds = 0, nanoseconds = 0, calibration = 0, prev_time = 0, next_time = 0;
  do{
    getline(run,line);
  }
  while (line.substr(11,5) != "EVENT");
  // for (int n = 0; n <= 108;n++) {getline(run,line);}

  //   for (k = 0; k <= 50000; k++){

   do  {k+= 1;// cout << "INIZIO DEL DO" << endl;
	 //  if (m%5000 == 0) cout << m << " eventi analizzati..." << endl;
     if (line.find("EVENT") > 15 ) {getline(run,line);continue;}//Durante il run compaiono righe non di evento, se non lo trova find restituisce un numero molto grande
    ch1 = 0;
    ch2 = 0;
    ch3 = 0;

    linecount = GetHitCount(line);

    if (linecount == 0){getline(run,line); continue;}

       point* hit = new point[linecount+3];//alloca la memoria per tutti i punti dell'evento
       //   cout << point::n << endl;
  // cout << GetHitCount(line) << endl;
  for (unsigned int i = 0; i < line.size()+1; i++) {
    if(!isalnum(line[i]) && line[i] != '.' && line[i] != '-') {
      j = j+1;

      if (entry != "" && j == 4) {stringstream(entry) >> number; cout << "EVENTO NUMERO: " << number << endl;} 
      if (entry != "" && j == 5) {stringstream(entry) >> seconds;}
      if (entry != "" && j == 6) {stringstream(entry) >> nanoseconds;}
      if (entry != "" && j == 7) {

	stringstream(entry) >> calibration;
	next_time = seconds*1000 + nanoseconds*4E7/calibration/1E6; //ms
	if (prev_time != 0) 
	  {time_dist->Fill(next_time-prev_time);} 
	
	//	cout << (next_time == prev_time) << " " << next_time << " " <<  prev_time << " " << calibration << endl;
	prev_time = next_time;
}
        
      if (entry != "" && j >= 9) {

	stringstream(entry) >> number;
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


  double* x = new double[linecount];
  double* y = new double[linecount];
  double* z = new double[linecount];

   for (int q = 0; q < linecount; q++){
     x[q] = hit[q].x;
     y[q] = hit[q].y;
     z[q] = hit[q].z;}


   TGraph2D* evdisplay = new TGraph2D(linecount,x,y,z);
    evdisplay->SetTitle("Event Display;X;Y;Z");

    // qui inserire event display
    //         if (k == 148) evdisplay->Write(); 

  //Plot distanze tra due hit, per le camere 2 e 3 aggiunto un offset di ch1 o (ch1+ch2) perché gli 
  // hit sono in ordine di camera nel file
  for (int h = ch1; h > 0; h--){
    for(int p = ch1-1; p < h && p >= 0; p--){
            
      d1 = (pow(pow(hit[h].x-hit[p].x,2)+pow(hit[h].y-hit[p].y,2),0.5)+d1*k1)/(k1+1); 
      disdist1->Fill(pow(pow(hit[h].x-hit[p].x,2)+pow(hit[h].y-hit[p].y,2),0.5)); 

      dx1 = (hit[h].x-hit[p].x+dx1*k1)/(k1+1); disdistx1->Fill(hit[h].x-hit[p].x);

      dy1 = (hit[h].y-hit[p].y+dy1*k1)/(k1+1); disdisty1->Fill(hit[h].y-hit[p].y);
      k1++;
    }}

  av_dist1->Fill(d1);
  av_distx1->Fill(dx1);
  av_disty1->Fill(dy1);

  for (int h = ch2; h > 0; h--){
    for(int p = ch2-1; p < h && p >= 0; p--){
            
      d2 = (pow(pow(hit[h+ch1].x-hit[p+ch1].x,2)+pow(hit[h+ch1].y-hit[p+ch1].y,2),0.5)+d2*k2)/(k2+1); 
      disdist2->Fill(pow(pow(hit[h+ch1].x-hit[p+ch1].x,2)+pow(hit[h+ch1].y-hit[p+ch1].y,2),0.5)); 

      dx2 = (hit[h+ch1].x-hit[p+ch1].x+dx2*k2)/(k2+1); disdistx2->Fill(hit[h+ch1].x-hit[p+ch1].x);

      dy2 = (hit[h+ch1].y-hit[p+ch1].y+dy2*k2)/(k2+1); disdisty2->Fill(hit[h+ch1].y-hit[p+ch1].y);
      k2++;
    }}

  av_dist2->Fill(d2);
  av_distx2->Fill(dx2);
  av_disty2->Fill(dy2);

  for (int h = ch3; h > 0; h--){
    for(int p = ch3-1; p < h && p >= 0; p--){
            
      d3 = (pow(pow(hit[h+ch1+ch2].x-hit[p+ch1+ch2].x,2)+pow(hit[h+ch1+ch2].y-hit[p+ch1+ch2].y,2),0.5)+d3*k3)/(k3+1); 
      disdist3->Fill(pow(pow(hit[h+ch1+ch2].x-hit[p+ch1+ch2].x,2)+pow(hit[h+ch1+ch2].y-hit[p+ch1+ch2].y,2),0.5)); 

      dx3 = (hit[h+ch1+ch2].x-hit[p+ch1+ch2].x+dx3*k3)/(k3+1); disdistx3->Fill(hit[h+ch1+ch2].x-hit[p+ch1+ch2].x);

      dy3 = (hit[h+ch1+ch2].y-hit[p+ch1+ch2].y+dy3*k3)/(k3+1); disdisty3->Fill(hit[h+ch1+ch2].y-hit[p+ch1+ch2].y);
      k3++;
    }}	

  av_dist3->Fill(d3);
  av_distx3->Fill(dx3);
  av_disty3->Fill(dy3);

  d1 = 0; d2 = 0; d3 = 0; dx1 = 0; dx2 = 0; dx3 = 0; dy1 = 0; dy2 = 0; dy3 = 0;
  k1 = 0; k2 = 0; k3 = 0;
  

  //   cout << "DOPO FOR" << endl;
 //Controllo la bontà dell'evento: ha almeno un hit per camera?
  if (ch1 < 1 || ch2 < 1 || ch3 < 1){/*cout << "EVENTO NON BUONO" << endl;*/ 
    j = 0; 
    delete [] hit;
    // cout << point::n << endl;
    delete[] x;
    delete[] y;
    delete[] z;
     delete evdisplay;
    getline(run,line);
    entry = ""; 

    continue;} //Evento non buono: prossimo evento

  else { //evento con almeno un hit per camera
   
    //    for (int kk = 0; kk < 3; kk++) {
    //	cout << hit[kk].x << endl;
    //	cout << hit[kk].y << endl;
    //	cout << hit[kk].z << endl;}
   
 manyhits2->Fill(ch2,ch1+ch2+ch3);

   for (int a = 0; a < (ch1); a++) {
     for (int b = 0; b < (ch2); b++) {
       for (int c = 0; c < (ch3); c++) {
   	 	 n1.SetPoints(hit[a],hit[b+ch1],hit[c+ch1+ch2]); //considero tutte le combinazioni di triplette
   		  //		  cout << "PRIMA FIT" << endl;
   		  n1.Fit();
		  if (n1.yvert || n1.zvert){ cout << "Verticalità a k: " << k << endl;
		    getline(cin,test);}
		  
		  //	  n1.XZFit();
   		  //  cout << "DOPO FIT" << endl;
		  //	  tempchi = (n1.YZGetChisquare()+n1.ZXGetChisquare()+n1.YXGetChisquare())/3;
		  tempchi = (n1.YXGetChisquare_m() + n1.YZGetChisquare_m() + n1.ZXGetChisquare_m())/3;
   		  //		  cout << "DOPO CHI" << endl;
   		  if (tempchi < chi && tempchi != 0 && !n1.yzfr) {chi = tempchi; yxtempchi = n1.YXGetChisquare_m(); yztempchi = n1.YZGetChisquare_m(); zxtempchi = n1.ZXGetChisquare_m(); theta = n1.GetTheta(); phi = n1.GetPhi();/*parameter = n1.YZGetParameter(1);besty = n1.yv;*/ besta = a; bestb = b; bestc = c; bestvalid = (n1.yxfr || n1.zxfr || n1.yzfr); bestxvert = n1.xvert;} else continue;
       }}}
   // if (parameter > 100000) {cout << "THETA: " << theta << endl; cout << "PHI: " << phi << endl; cout << "y SOSPETTE: " << besty[0] << " " << besty[1] << " " << besty[2] << endl; cout << "k: " << k << endl;}
   if (bestxvert) {cout << "VERTICALE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl; cout << phi << endl; cout << "valido: " << bestvalid << endl;  cout << "xyChi: " << yxtempchi << endl; cout << "yzchi: " << yztempchi << endl; cout << "zxchi: " << zxtempchi << endl;}


   disxy1->Fill(n1.GetCoordinate(0,0),n1.GetCoordinate(1,0));
   disxy2->Fill(n1.GetCoordinate(0,1),n1.GetCoordinate(1,1));
   disxy3->Fill(n1.GetCoordinate(0,2),n1.GetCoordinate(1,2));
   
   disx1->Fill(n1.GetCoordinate(0,0));
   disx2->Fill(n1.GetCoordinate(0,1));
   disx3->Fill(n1.GetCoordinate(0,2));

   disy1->Fill(n1.GetCoordinate(1,0));
   disy2->Fill(n1.GetCoordinate(1,1));
   disy3->Fill(n1.GetCoordinate(1,2));

   
  
   disyxchi->Fill(yxtempchi);
   diszxchi->Fill(zxtempchi);
   disyzchi->Fill(yztempchi);
      dischi->Fill(chi);
  
   //Riempiamo gli istogrammi di theta e phi se il fit è andato bene. Se è verticale considero solamente una sezione
      //     if (phi > 1.5 && phi < 1.64 && yztempchi > 30 && bestvert) {cout << "y sospette: " << besty[0] << '\t' << besty[1] << '\t' << besty[2] << endl; cout << "k: " << k << endl; cin.get();}
      if( (yxtempchi < 6.63) && (yztempchi < 6.63) && (zxtempchi < 6.63) && !bestvalid){
     // cout << "Fit con chi2: " << chi << endl;
     // cout << "Theta: " << n1.GetTheta() << endl;
     // cout << "Phi: " << n1.GetPhi() << endl;
	//   if (theta <= 0.025){cout << "TROVATO THETA = 0 in evento" << k << " con theta: " << theta << " e phi: " << phi << " xyparamter: " << n1.XYGetParameter(1) << " yzparameter: " << n1.YZGetParameter(1) << endl; cin.get();}   
	//	if (theta < 1E-4) {cout << "y SOSPETTE: " << besty[0] << " " << besty[1] << " " << besty[2] << endl; cout << "theta acc: " << theta << endl; cout << "k da contr: " << k << endl; cin.get();}
	//	if (theta > 1.22) {cout << "HUGE THETA at k: " << k << endl;}
	if (theta*180/3.14159 <= 0.1) {piccolitheta->Fill(phi*180/3.141592654,theta*180/3.141592654); cout << "Verticalità: " << bestxvert << endl; cout << theta*180/3.141592654 << endl; cout << phi*180/3.141592654 << endl; cout << "besta,b,c: " << besta << " " << bestb << " " << bestc << endl;cout << "CHI " << chi << endl;}
	disprimatheta->Fill(theta*180/3.141592654);
	disprimaphi->Fill(phi*180/3.141592654);
	primathetaphi->Fill(phi*180/3.141592654,theta*180/3.141592654);
      }
 


      hit[linecount].SetValues(hit[besta].x,hit[besta].y,hit[besta].z);
      hit[linecount+1].SetValues(hit[bestb+ch1].x,hit[bestb+ch1].y,hit[bestb+ch1].z);
      hit[linecount+2].SetValues(hit[bestc+ch1+ch2].x,hit[bestc+ch1+ch2].y,hit[bestc+ch1+ch2].z);

      if(hit[linecount].z == hit[linecount+1].z || hit[linecount+1].z == hit[linecount+2].z) {cout << "ERROREEEEEEEEEEEEEEEEEEE a k = " << k << endl; cout << "CHI: " << chi << endl; cout << "THETA,PHI: " << theta*180/3.14159 << " " << phi*180/3.14159 << endl; getline(cin,test);}
      
      // n1.SetPoints(hit[linecount],hit[linecount+1],hit[linecount+2]); //tripletta punti migliori
      // n1.XYFit();
      // n1.YZFit();
      
      // XYint = n1.XYGetInter


  }
  j = 0;
    delete[] x;
    delete[] y;
    delete[] z;
     delete evdisplay;
  delete[] hit;
  // cout << point::n << endl;
  getline(run,line);
  entry = "";
  chi = 1E30;
  tempchi = 0;
  //Bin per chamber histograms
   hpc1->Fill(ch1);
   hpc2->Fill(ch2);
   hpc3->Fill(ch3);
   disz->SetBinContent(4,disz->GetBinContent(4)+ch1);
   disz->SetBinContent(8,disz->GetBinContent(8)+ch2);
   disz->SetBinContent(12,disz->GetBinContent(12)+ch3);
   
   // cout << "FINE DEL DO" << endl;
	  }



   while (!run.eof());

   double yxchi_integral = 0;
   for (int i = 81; i <= 300; i++) {
     yxchi_integral += disyxchi->GetBinContent(i);
   }

   cout << "Integrale da 8 a 10 di yxchi: " << yxchi_integral << endl;
   cout << "Dovrebbero essere: " << disyxchi->GetEntries()*0.0046 << endl;
   

     // TCanvas* thetacanv = new TCanvas();
     // thetacanv->SetGrid();
     // thetacanv->cd();
     // distheta->Draw();
  
      //TCanvas* phicanv = new TCanvas();
       //phicanv->SetGrid();
      //phicanv->cd();
  //   disphi->Draw();
   manyhitsline->Write();
   manyhits2->Write();
   primathetaphi->Write();
   disdist1->Write();
   disdist2->Write();
   disdist3->Write();
   disdistx1->Write();
   disdistx2->Write();
   disdistx3->Write();
   disdisty1->Write();
   disdisty2->Write();
   disdisty3->Write();
   disxy1->Write();
   disxy2->Write();
   disxy3->Write();
   disx1->Write();
   disx2->Write();
   disx3->Write();
   disy1->Write();
   disy2->Write();
   disy3->Write();
   disz->Write();
    hpc1->Write();
    hpc2->Write();
    hpc3->Write();
    dischi->Write();
    disyxchi->Write();
    diszxchi->Write();
    disyzchi->Write();
    av_dist1->Write();
    av_dist2->Write();
    av_dist3->Write();
    av_distx1->Write();
    av_distx2->Write();
    av_distx3->Write();
    av_disty1->Write();
    av_disty2->Write();
    av_disty3->Write();
    
     disprimatheta->Write();
     disprimaphi->Write();
     piccolitheta->Write();
     time_dist->Write();
     rfile.Close();
     run.close();
}
  
