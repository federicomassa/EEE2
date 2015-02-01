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
  for (unsigned int i = 0; i < str.size()+1; i++) {
    if(!isalnum(str[i]) && str[i] != '.' && str[i] != '-') //se trova un carattere non alfanumerico o . allora è la fine di un blocco
      s = s+1;
  }
  return (s-9)/3;
}; // 9 blocchi prima dell'inizio delle 3 coordinate



void tracks(){  
  bool bestvert = false;
  bool clustvert = false;
  //  double *besty = new double[3];
  // ifstream run("../Data/EEE_Prova_topbottom8900__20140530_174533.txt"); //INPUT FILE

  ifstream run("../EEEData/CORR_EEE_PISA01TestRun4Telescopes_20140507_014455.txt"); //INPUT FILE
  //  int point::n = 0;
    // TFile rfile("Disttopbot.root","RECREATE");
 TFile rfile("Dist_all_trasl.root","RECREATE");
 TF1* manyhitsline = new TF1("prova","[0]+[1]*x",0,15);
  TH2F* manyhits2 = new TH2F("manyhits","Many-hits events correlation: Ch 2", 30,0,15,100,0,50);
  TH1F* dischi = new TH1F("dischi","Chi2 distribution; chi2; #", 100,0,10);
  TH1F* hpc1 = new TH1F("hpc1", "Hit per chamber / Chamber 1; #Hits;# ", 20,0,20);
  TH1F* hpc2 = new TH1F("hpc2", "Hit per chamber / Chamber 2;#Hits;#", 20,0,20);
  TH1F* hpc3 = new TH1F("hpc3", "Hit per chamber /Chamber 3;#Hits;#", 20,0,20);
  TH1F* distheta = new TH1F("dist","Theta distribution;Theta(deg);#", 50, 0,90);
  TH1F* disphi = new TH1F("disp","Phi distribution;Phi(deg);#", 90, 0, 360);
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
  TH2F* thetaphi = new TH2F("phi-teta","Phi-Theta Correlation;Phi(deg);Theta(deg)", 90,0,360,50,0,90);
  TH2F* preposttheta = new TH2F("prepostth","Theta PRE Cluster vs POST Cluster;Theta PRE (deg);Theta POST (deg)", 90,0,90,90,0,90);
  TH2F* prepostphi = new TH2F("prepostphi","Phi PRE Cluster vs POST Cluster;Phi PRE (deg);Phi POST (deg)", 180,0,360,180,0,360);
  TH1F* disprimatheta = new TH1F("disprimatheta","Distribuzione theta prima del cluster; Theta (gradi); #",90,0,90);
  TH1F* disprimaphi = new TH1F("disprimaphi","Distribuzione phi prima del cluster; Phi (gradi);#",180,0,360);
  triplet n1;
  double chi = 1000, tempchi = 0,xytempchi = 0, /*xztempchi = 0, */ yztempchi = 0, theta = 1000, phi = 1000;
  int j = 0;
  double number = 0;
  double XYint,YZint,XYpen,YZpen,XZpen,midx1,midx2,midx3,midy1,midy2,midy3, xav1=0, xav2=0, xav3=0, yav1=0, yav2=0, yav3=0;
    double jx1=0;
    double jx2=0;
    double jy1=0;
    double jy2=0;
    double jx3=0;
    double jy3=0;
    double xyrad=0;
    int besta=100, bestb=100,bestc=100;
    double primatheta,primaphi;



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

  //   for (k = 0; k <= 50000; k++){

   do  {k+= 1;// cout << "INIZIO DEL DO" << endl;
	 //  if (m%5000 == 0) cout << m << " eventi analizzati..." << endl;
     if (line.find("EVENT") > 15 ) {getline(run,line);continue;}//Durante il run compaiono righe non di evento, se non lo trova find restituisce un numero molto grande
    ch1 = 0;
    ch2 = 0;
    ch3 = 0;

    linecount = GetHitCount(line);
    //    cout << line << endl;
    if (linecount == 0){getline(run,line); continue;}

       point* hit = new point[linecount+3];//alloca la memoria per tutti i punti dell'evento
       //   cout << point::n << endl;
  // cout << GetHitCount(line) << endl;
  for (unsigned int i = 0; i < line.size()+1; i++) {
    if(!isalnum(line[i]) && line[i] != '.' && line[i] != '-') {
      j = j+1;
      if (entry != "" && j == 4) {stringstream(entry) >> number; cout << "EVENTO NUMERO: " << number << endl;} 
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
	  disdist1->Fill(pow(pow(hit[h].x-hit[p].x,2)+pow(hit[h].y-hit[p].y,2),0.5));
    disdistx1->Fill(hit[h].x-hit[p].x);
    disdisty1->Fill(hit[h].y-hit[p].y);
    }}

 for (int h = ch2; h > 0; h--){
    for(int p = ch2-1; p < h && p >= 0; p--){
	  disdist2->Fill(pow(pow(hit[h+ch1].x-hit[p+ch1].x,2)+pow(hit[h+ch1].y-hit[p+ch1].y,2),0.5));
	      disdistx2->Fill(hit[h+ch1].x-hit[p+ch1].x);
    disdisty2->Fill(hit[h+ch1].y-hit[p+ch1].y);
    }}

 for (int h = ch3; h > 0; h--){
    for(int p = ch3-1; p < h && p >= 0; p--){
	  disdist3->Fill(pow(pow(hit[h+ch1+ch2].x-hit[p+ch1+ch2].x,2)+pow(hit[h+ch1+ch2].y-hit[p+ch1+ch2].y,2),0.5));
    disdistx3->Fill(hit[h+ch1+ch2].x-hit[p+ch1+ch2].x);
    disdisty3->Fill(hit[h+ch1+ch2].y-hit[p+ch1+ch2].y);
    }}
	  

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
   		  n1.XYFit();
   		  n1.YZFit();
		  //	  n1.XZFit();
   		  //  cout << "DOPO FIT" << endl;
		  if (!n1.vert)
		  tempchi = (n1.XYGetChisquare_m()/*+n1.XZGetChisquare_m()*/ +n1.YZGetChisquare_m())/2; //tolto xz per sicurezza: potrebbe essere verticale
		  else tempchi = n1.YZGetChisquare()/2;
   		  //		  cout << "DOPO CHI" << endl;
   		  if (tempchi <= chi && tempchi != 0) {chi = tempchi; xytempchi = n1.XYGetChisquare_m(); yztempchi = n1.YZGetChisquare_m(); /* xztempchi = n1.XZGetChisquare_m();*/ theta = n1.GetTheta(); primatheta = theta; disprimatheta->Fill(theta*180/3.14159); phi = n1.GetPhi(); primaphi = phi;disprimaphi->Fill(phi*180/3.14159);/*parameter = n1.YZGetParameter(1);besty = n1.yv;*/ bestvert = n1.vert; besta=a;bestb=b;bestc=c;}
       }}}
   // if (parameter > 100000) {cout << "THETA: " << theta << endl; cout << "PHI: " << phi << endl; cout << "y SOSPETTE: " << besty[0] << " " << besty[1] << " " << besty[2] << endl; cout << "k: " << k << endl;}
   disxy1->Fill(n1.GetCoordinate(0,0),n1.GetCoordinate(1,0));
   disxy2->Fill(n1.GetCoordinate(0,1),n1.GetCoordinate(1,1));
   disxy3->Fill(n1.GetCoordinate(0,2),n1.GetCoordinate(1,2));
   
   disx1->Fill(n1.GetCoordinate(0,0));
   disx2->Fill(n1.GetCoordinate(0,1));
   disx3->Fill(n1.GetCoordinate(0,2));

   disy1->Fill(n1.GetCoordinate(1,0));
   disy2->Fill(n1.GetCoordinate(1,1));
   disy3->Fill(n1.GetCoordinate(1,2));
   
      dischi->Fill(chi);

      {//C                         l                            u                        s                     t                  e
	//ricostruiamo la tripletta dal punto migliore

	n1.SetPoints(hit[besta],hit[ch1+bestb],hit[ch1+ch2+bestc]);
	cout << "z1: " << hit[besta].z << " z2: " << hit[ch1+bestb].z << " z3: " << hit[ch1+ch2+bestc].z << endl;
	cout << "besta: " << besta << " bestb" << bestb << " bestc " << bestc << endl;
	

	XYint = n1.XYGetParameter(0);
	 YZint = n1.YZGetParameter(0);
	 //	 XZint = n1.XZGetParameter(0);
 
	 XYpen= n1.XYGetParameter(1);
	 //	 cout << XYpen << endl;
	 YZpen= n1.YZGetParameter(1);
	 //XZpen = n1.XZGetParameter(1);

	 midx1 = (zch1-YZint-XYint*YZpen)/(XYpen*YZpen);
	 midx2 = (zch2-YZint-XYint*YZpen)/(XYpen*YZpen);
	 midx3 = (zch3-YZint-XYint*YZpen)/(XYpen*YZpen);

	 midy1 = (zch1-YZint)/YZpen;
	 midy2 = (zch2-YZint)/YZpen;
	 midy3 = (zch3-YZint)/YZpen;

   // 	 for(int m1 = 0; m1 < ch1; m1++) {
   // 	   if( (hit[m1].x-midx1)*(hit[m1].x-midx1)+(hit[m1].y-midy1)*(hit[m1].y-midy1)  <= xyrad*xyrad) {
   // 	     xav1 = (xav1*jx1+ hit[m1].x)/(jx1+1);
   // 	     yav1 = (yav1*jy1+ hit[m1].y)/(jy1+1);jy1+= 1;
   // 	     jx1+= 1;
   // 	   }
	  
   // 	 }
   // 	 if (jx1==0){xav1 = hit[besta].x; yav1=hit[besta].y;}
	  
	 
   // //if(absval(hit[m1].y - midy1) <= yrad1) {yav1 = (yav1*jy1+ hit[m1].y)/(jy1+1);jy1+= 1;}
   // 	     // else {yav1 = hit[besta].y;}
	  
   // 	 for(int m2 = ch1; m2 < ch1+ch2; m2++) {
   // 	   if( (hit[m2].x-midx2)*(hit[m2].x-midx2)+(hit[m2].y-midy2)*(hit[m2].y-midy2)  <= xyrad*xyrad) {
   // 	     cout << "primissima" << endl;
   // 	     xav2 = (xav2*jx2+ hit[m2].x)/(jx2+1);
   // 	     yav2 = (yav2*jy2+ hit[m2].y)/(jy2+1);jy2+= 1;
   // 	     jx2+= 1;

   // 	   }

   // 	 }
   // 	 if(jx2==0) {xav2 = hit[bestb+ch1].x;yav2=hit[bestb+ch1].y;}

   // 	 for(int m3 = ch2+ch1; m3 < ch3+ch1+ch2; m3++) {
   // 	   if( (hit[m3].x-midx3)*(hit[m3].x-midx3)+(hit[m3].y-midy3)*(hit[m3].y-midy3)  <= xyrad*xyrad) {
   // 	     xav3 = (xav3*jx3+ hit[m3].x)/(jx3+1);
   // 	     yav3 = (yav3*jy3+ hit[m3].y)/(jy3+1);jy3+= 1;
   // 	     jx3+= 1;
   // 	   }

   // 	 }
   // 	 if (jx3==0) {xav3 = hit[bestc+ch1+ch2].x;yav3=hit[bestc+ch1+ch2].y;}	
	 
//   
	 jx1 = 0;
	 jy1 = 0;
	 jx2 = 0;
	 jy2 = 0;
	 jx3 = 0;
	 jy3 = 0;

	 // for(int m2 = 0; m2 < ch2; m2++) {
	 //   if(absval(hit[m2+ch1].x - midx2) <= xrad2) {xav2 = (xav2*jx2+ hit[m2+ch1].x)/(jx2+1);jx2+= 1;}
	 //   else {xav2 = hit[bestb+ch1].x;}
	   
	 //    if(absval(hit[m2+ch1].y - midy2) <= yrad2) {yav2 = (yav2*jy2+ hit[m2+ch1].y)/(jy2+1);jy2+= 1;}
	 //   else {yav2 = hit[bestb+ch1].y;}
	 // }

	 // for(int m3 = 0; m3 < ch3; m3++) {
	 //   if(absval(hit[m3+ch1+ch2].x - midx3) <= xrad3) {xav3 = (xav3*jx3+ hit[m3+ch1+ch2].x)/(jx3+1);jx3+= 1;}
	 //   else {xav3 = hit[bestc+ch1+ch2].x;}
	   
	 //    if(absval(hit[m3+ch1+ch2].y - midy3) <= yrad3) {yav3 = (yav3*jy3+ hit[m3+ch1+ch2].y)/(jy3+1);jy3+= 1;}
	 //   else {yav3 = hit[bestc+ch1+ch2].y;}
	 // }


	 hit[linecount].SetValues(hit[besta].x,hit[besta].y,hit[besta].z);
	 hit[linecount+1].SetValues(hit[bestb+ch1].x,hit[bestb+ch1].y,hit[bestb+ch1].z);
	 hit[linecount+2].SetValues(hit[ch1+ch2+bestc].x,hit[ch1+ch2+bestc].y,hit[ch1+ch2+bestc].z);

  n1.SetPoints(hit[linecount],hit[linecount+1],hit[linecount+2]);
  n1.XYFit();n1.YZFit();
  chi = (n1.XYGetChisquare_m()+n1.YZGetChisquare_m())/2;
      } //END CLUSTER
  
      clustvert = n1.vert;

   //Riempiamo gli istogrammi di theta e phi se il fit è andato bene. Se è verticale considero solamente una sezione
      //     if (phi > 1.5 && phi < 1.64 && yztempchi > 30 && bestvert) {cout << "y sospette: " << besty[0] << '\t' << besty[1] << '\t' << besty[2] << endl; cout << "k: " << k << endl; cin.get();}
      if( ((xytempchi>0 || clustvert) && yztempchi/**xztempchi*/ > 0)/* && (chi<5)*/  && (xytempchi < 10 || clustvert) && (yztempchi < 10)/* && (xztempchi < 10)*/ ){
     // cout << "Fit con chi2: " << chi << endl;
      if(clustvert)      cout << "VERTICALEEEEEEEEEEEEEEE" << endl;
	cout << "Theta: " << n1.GetTheta() << " Prima era " << primatheta <<  endl;
     // cout << "Phi: " << n1.GetPhi() << endl;
	//   if (theta <= 0.025){cout << "TROVATO THETA = 0 in evento" << k << " con theta: " << theta << " e phi: " << phi << " xyparamter: " << n1.XYGetParameter(1) << " yzparameter: " << n1.YZGetParameter(1) << endl; cin.get();}   
	//	if (theta < 1E-4) {cout << "y SOSPETTE: " << besty[0] << " " << besty[1] << " " << besty[2] << endl; cout << "theta acc: " << theta << endl; cout << "k da contr: " << k << endl; cin.get();}
	//	if (theta > 1.22) {cout << "HUGE THETA at k: " << k << endl;}
	theta = n1.GetTheta();
	phi = n1.GetPhi();
	preposttheta->Fill(primatheta*180/3.14159,theta*180/3.1415);
	prepostphi->Fill(primaphi*180/3.14159,phi*180/3.1415);
	thetaphi->Fill(phi*180/3.14159,theta*180/3.14159);
	distheta->Fill(theta*180/3.14159);
	disphi->Fill(phi*180/3.14159);
	if (theta<0.0001){cout << "chi " << chi << " midx1,2,3 - midy1,2,3 " << midx1 << " " << midx2 << " " << midx3 << " " << midy1 << " " << midy2 << " " << midy3 << endl;}
       }
   //  else {
   //  jj += 1;
   //  if (jj == 1) {chi2 << chi << endl; n1.XYDraw(); n1.YZDraw();}
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
  chi = 1000;
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



       	  while (k <= 10000);
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
   thetaphi->Write();
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
   preposttheta->Write();
   prepostphi->Write();
    hpc1->Write();
    hpc2->Write();
    hpc3->Write();
    dischi->Write();
     distheta->Write();
     disphi->Write();
     disprimatheta->Write();
     disprimaphi->Write();
     rfile.Close();
     run.close();
}
  
