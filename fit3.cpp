//rispetto a fit.cpp cambio gli assi in modo da non avere problemi quando
// x = cost

#include "hit3.cpp"
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLinearFitter.h>

using namespace std;

double absval(double d) {
  if (d >= 0) return (d);
  else return (-d);}

double sqrt(double a) {return TMath::Sqrt(a);}
double sqr(double a) {return (a*a);}

class triplet{
public:

  //calcolo "manuale" del chi2 per i fit che danno particolari problemi (picco in theta = 0..)
  double YXGetChisquare_m() {
    YXFit();
    
    double eqerr[3];
    double result = 0;
    // somma in quadratura: attenzione: viene usata la pendenza tra due punti per il calcolo di eqerr
        for (unsigned n = 0; n < 3; n++) eqerr[n] = sqrt(sqr(xerr[n]) + sqr(YXGetParameter(1)*yerr[n]));

    for (unsigned n = 0; n < 3; n++) result += sqr((xv[n]- (YXGetParameter(0) + YXGetParameter(1)*yv[n]))/eqerr[n]);
    return result; 

  }

  double ZXGetChisquare_m() {
    ZXFit();
    double eqerr[3];
    double result = 0;
    // somma in quadratura: attenzione: viene usato XYGetParameter(1) per il calcolo di eqerr, possibile pericolo
    for (unsigned n = 0; n < 3; n++) eqerr[n] = sqrt(sqr(xerr[n]) + sqr(ZXGetParameter(1)*zerr[n]));
    for (unsigned n = 0; n < 3; n++) result += sqr((xv[n]- (ZXGetParameter(0) + ZXGetParameter(1)*zv[n]))/eqerr[n]);
    return result;
  }

  double YZGetChisquare_m() {
    YZFit();
    double eqerr[3];
    double result = 0;
    // somma in quadratura
     for (unsigned n = 0; n < 3; n++) eqerr[n] = sqrt(sqr(zerr[n]) + sqr(YZGetParameter(1)*yerr[n]));
     //  for (unsigned n = 0; n < 3; n++) eqerr[n] = zerr[n] + YZGetSlope()*yerr[n];
    for (unsigned n = 0; n < 3; n++) result += sqr((zv[n]- (YZGetParameter(0) + YZGetParameter(1)*yv[n]))/eqerr[n]);
    return result;
  }


  //Calcolo intercetta e pendenza sui primi due punti per migliore inizializzazione dei parametri del fit
 double YXGetIntercept() {
   return (xv[0]-(xv[1]-xv[0])/(yv[1]-yv[0])*yv[0]) ;}
 double ZXGetIntercept() {
   return (xv[0]-(xv[1]-xv[0])/(zv[1]-zv[0])*zv[0]) ;}
 double YZGetIntercept() {
   return (zv[0]-(zv[1]-zv[0])/(yv[1]-yv[0])*yv[0]) ;}

 double YXGetSlope() {
    return (xv[1]-xv[0])/(yv[1]-yv[0]) ;}
 double ZXGetSlope() {
    return (xv[1]-xv[0])/(zv[1]-zv[0]) ;}
 double YZGetSlope() {
    return (zv[1]-zv[0])/(yv[1]-yv[0]) ;}

  double xv[3];
  double yv[3];
  double zv[3];
  double xerr[3];
  double yerr[3];
  double zerr[3];
  Int_t yxfr;
  TF1 *yxfitfunc;
  TGraphErrors* yxgraph;
  Int_t zxfr;
  TF1* zxfitfunc;
  TGraphErrors* zxgraph;
  Int_t yzfr;
  TF1* yzfitfunc;
  TGraphErrors* yzgraph;
  double theta,phi;
  bool xvert, yvert,zvert;
  //public:
  triplet(){
    yxfr = 0;
    yzfr = 0;
    zxfr = 0;
    yxfitfunc = 0;
    zxfitfunc = 0;
    yzfitfunc = 0;
    yxgraph = 0;
    zxgraph = 0;
    yzgraph = 0;
    for (int i = 0; i < 3; i++) xv[i] = 0;
    for (int i = 0; i < 3; i++) yv[i] = 0;
    for (int i = 0; i < 3; i++) zv[i] = 0;
    theta = 0;
    phi = 0;
    for (int i = 0; i < 3; i++) xerr[i] = 0;
    for (int i = 0; i < 3; i++) yerr[i] = 0;
    for (int i = 0; i < 3; i++) zerr[i] = 0;
  };

  double YXGetChisquare(){return (yxfitfunc->GetChisquare());};
  double YXGetParameter(int i){return (yxfitfunc->GetParameter(i));};
  double ZXGetChisquare(){return (zxfitfunc->GetChisquare());};
  double ZXGetParameter(int i){return (zxfitfunc->GetParameter(i));};
  double YZGetChisquare(){return (yzfitfunc->GetChisquare());};
  double YZGetParameter(int i){return (yzfitfunc->GetParameter(i));};
  double GetPhi(){
    if (!xvert){
    phi = atan(1/yxfitfunc->GetParameter(1));
    if (yzfitfunc->GetParameter(1) > 0 && phi < 0) phi = phi + 4*atan(1.);
    if (yzfitfunc->GetParameter(1) < 0 && phi > 0) phi = phi + 4*atan(1.);
    if (yzfitfunc->GetParameter(1) < 0 && phi < 0) phi = phi + 8*atan(1.);
//Uso l'altra sezione per risolvere l'indecisione di 180° su phi
    return phi;}
    else {
      if (YZGetParameter(1) > 0) return 2*atan(1);
      else {return 6*atan(1);} 
    
    }
  }
  
 
  
  double GetTheta(){
    phi = GetPhi();
    theta = absval(atan(1/(yzfitfunc->GetParameter(1)*(sin(phi)))));
    return theta;};

void SetPoints(point x1, point x2, point x3){
    xv[0] = x1.x; xv[1] = x2.x; xv[2] = x3.x;
    yv[0] = x1.y ; yv[1] = x2.y; yv[2] = x3.y;
    zv[0] = x1.z; zv[1] = x2.z; zv[2] = x3.z;
    //controllo di verticalità
    if (xv[0] == xv[1] && xv[0] == xv[2]) xvert = true; else xvert = false;
    if (yv[0] == yv[1] && yv[0] == yv[2]) yvert = true; else yvert = false;
    if (zv[0] == zv[1] && zv[0] == zv[2]) zvert = true; else zvert = false;
      for (int i = 0; i < 3;i++) {xerr[i] = 0.8; yerr[i] = 1.5;zerr[i] = 0.5;}
      // for (int i = 0; i < 3;i++) {xerr[i] = 1; yerr[i] = 0;zerr[i] = 3.333;}
};  //incertezze di default, rivedere perché 2.
  
   triplet(point x1, point x2, point x3){
     xv[0] = x1.x; xv[1] = x2.x; xv[2] = x3.x;
     yv[0] = x1.y ; yv[1] = x2.y; yv[2] = x3.y;
     zv[0] = x1.z; zv[1] = x2.z; zv[2] = x3.z;
     for (int i = 0; i < 3;i++) {xerr[i] = 1.44; yerr[i] = 2;zerr[i] = 0.5;}};  //incertezze di default

     ~triplet() {
    delete yxfitfunc;
     delete zxfitfunc;
     delete yzfitfunc;
     delete yxgraph;
     delete zxgraph;
     delete yzgraph;};
    

  //  void SetErr(double ex[3],double ey[3], double ez[3]){xerr = ex; yerr = ey;zerr = ez;};

  void YXFit(){
   
    
    if (xerr[0] == 0 && xerr[1] == 0 && xerr[2] == 0 && yerr[0] == 0 && yerr[1] == 0 && yerr[2] == 0) {cout << "ERRORE: Incertezze non specificate" << endl;}
    else {
      //  TGraphErrors *graph1 = new TGraphErrors(3,x,y,xerr,yerr);
    //  TF1* fitfunc1 = new TF1("fittingfunction","[0]+[1]*x",-100,100);
    //   xyfitfunc = fitfunc1;
    // xygraph = graph1;
      delete yxfitfunc;//Dealloca prima di riallocarne uno nuovo
      delete yxgraph;

         double par0, par1;
        yxgraph = new TGraphErrors(3,yv,xv,yerr,xerr);
	yxgraph->LeastSquareLinearFit(3,par0,par1,yxfr,-400,400);
	yxfitfunc = new TF1("yxfitfunc","[0]+[1]*x",-400,400);
	if(yxfr == 0)  	yxfitfunc->SetParameters(par0, par1); 
	else yxfitfunc->SetParameters(1E10, 0);
	

        // yxgraph = new TGraphErrors(3,yv,xv,yerr,xerr); 
	// yxfr = yxgraph->Fit("1 ++ x","0QS");
	// yxfitfunc = yxgraph->GetFunction("1 ++ x");
	
    } 

  };
  
  // void XYDraw(){
  //   xygraph->SetLineColor(kBlue);
  //    xyfitfunc->SetLineColor(kRed);
  // TCanvas* XYlinearfit = new TCanvas();
  // XYlinearfit->SetTitle("XY Linear Fit");
  // XYlinearfit->SetGrid();
  // xygraph->GetXaxis()->SetTitle("x [cm]");
  // xygraph->GetYaxis()->SetTitle("y [cm]");
  // xygraph->GetXaxis()->CenterTitle();
  // xygraph->GetYaxis()->CenterTitle();
  // xygraph->DrawClone("APE");
  // xyfitfunc->DrawClone("SAME");
  // };

  void ZXFit(){   
    if (xerr[0] == 0 && xerr[1] == 0 && xerr[2] == 0 && zerr[0] == 0 && zerr[1] == 0 && zerr[2] == 0) {cout << "ERRORE: Incertezze non specificate" << endl;}
    else {
      //  TGraphErrors *graph1 = new TGraphErrors(3,x,y,xerr,yerr);
    //  TF1* fitfunc1 = new TF1("fittingfunction","[0]+[1]*x",-100,100);
    //   xyfitfunc = fitfunc1;
    // xygraph = graph1;
      delete zxfitfunc;  //Dealloca prima di riallocarne uno nuovo
      delete zxgraph;
     
        double par0, par1;

        zxgraph = new TGraphErrors(3,zv,xv,zerr,xerr);
	zxgraph->LeastSquareLinearFit(3,par0,par1,zxfr,-400,400);
	zxfitfunc = new TF1("zxfitfunc","[0]+[1]*x",-400,400);
	if(zxfr == 0)  	zxfitfunc->SetParameters(par0, par1); 
	else zxfitfunc->SetParameters(1E10, 0);

        // zxgraph = new TGraphErrors(3,zv,xv,zerr,xerr); 
	// zxfr = zxgraph->Fit("1 ++ x","0QS");
	// zxfitfunc = zxgraph->GetFunction("1 ++ x");
	
    }
  };
  
  // void XZDraw(){
  // xzgraph->SetLineColor(kBlue);
  // xzfitfunc->SetLineColor(kRed);
  // TCanvas* XZlinearfit = new TCanvas();
  // XZlinearfit->SetTitle("XZ Linear Fit");
  // XZlinearfit->SetGrid();
  // xzgraph->GetXaxis()->SetTitle("x [cm]");
  // xzgraph->GetYaxis()->SetTitle("z [cm]");
  // xzgraph->GetXaxis()->CenterTitle();
  // xzgraph->GetYaxis()->CenterTitle();
  // xzgraph->DrawClone("APE");
  // xzfitfunc->DrawClone("SAME");
  // };

void YZFit(){   
    if (yerr[0] == 0 && yerr[1] == 0 && yerr[2] == 0 && zerr[0] == 0 && zerr[1] == 0 && zerr[2] == 0) {cout << "ERRORE: Incertezze non specificate" << endl;}
    else {
      delete yzfitfunc; 
      delete yzgraph;

        double par0, par1;
        yzgraph = new TGraphErrors(3,yv,zv,yerr,zerr);
	yzgraph->LeastSquareLinearFit(3,par0,par1,yzfr,-400,400);
	yzfitfunc = new TF1("yzfitfunc","[0]+[1]*x",-400,400);
	if(yzfr == 0)  	yzfitfunc->SetParameters(par0, par1); 
	else yzfitfunc->SetParameters(1E10, 0);

    
        // yzgraph = new TGraphErrors(3,yv,zv,yerr,zerr); 

	// yzfitfunc = yzgraph->GetFunction("1 ++ x");
	
    
    }}
  
  // void YZDraw(){
  //   yzgraph->SetLineColor(kBlue);
  //    yzfitfunc->SetLineColor(kRed);
  // TCanvas* YZlinearfit = new TCanvas();
  // YZlinearfit->SetTitle("YZ Linear Fit");
  // YZlinearfit->SetGrid();
  // yzgraph->GetXaxis()->SetTitle("y [cm]");
  // yzgraph->GetYaxis()->SetTitle("z [cm]");
  // yzgraph->GetXaxis()->CenterTitle();
  // yzgraph->GetYaxis()->CenterTitle();
  // yzgraph->DrawClone("APE");
  // yzfitfunc->DrawClone("SAME");
  // };

  void Fit(){
    ZXFit();
    YXFit();
    YZFit();};

  double GetCoordinate(int coord, int pn) {
    switch (coord){
    case(0):
      return (xv[pn]);
      break;
    case(1):
      return (yv[pn]);
      break;
    case(2):
      return (zv[pn]);
      break;
    default:
      cout << "INVALID COORDINATE NUMBER: " << coord << endl;
      return (-5000);
    }};
    
    

};
    

void fit(){
  point* p1 = new point; p1->SetValues(-39.2917, 61.1502, 53.2); 
  point* p2 = new point; p2->SetValues(-39.2917, 63.5019, 0);
  point* p3 = new point; p3->SetValues(-35.875, 57.4314, -52.8);
  // point* p1 = new point; p1->SetValues(-30.00,138.24,145.00); 
  // point* p2 = new point; p2->SetValues(-32.00,148.55,85.00);
  // point* p3 = new point; p3->SetValues(-34.00,158.24,24.00);
  triplet n1; n1.SetPoints(*p1,*p2,*p3);
  n1.YXFit();
  n1.YZFit();
  n1.ZXFit();
  n1.yzgraph->DrawClone("APE");
  n1.yzfitfunc->DrawClone("same");
  cout << "Valido? " << n1.yzfr << endl;
  cout << "xy: " << n1.YXGetChisquare() << endl;
  cout << "xz: " << n1.ZXGetChisquare() << endl;
  cout << "yz: " << n1.YZGetChisquare() << endl;
  
  // TGraphErrors* g1 = new TGraphErrors(3,n1.yv,n1.zv,n1.yerr,n1.zerr);
  // cout << "intercetta: " << n1.YZGetIntercept();
  // cout << "pendenza: " << n1.YZGetSlope();
  // TF1* f1 = new TF1("lin","[0]+[1]*x",-100,100);
  // f1->SetParameters(n1.YZGetIntercept(),n1.YZGetSlope());
  // g1->Fit(f1,"0");
  // g1->Draw("APE");
  // f1->DrawCopy("same");
  // cout << "XY 0: " << n1.XYGetParameter(0) << endl;
  // cout << "XY 1: " << n1.XYGetParameter(1) << endl;
  // cout << "THETA: " << n1.GetTheta() << endl;
  // cout << "PHI: " << n1.GetPhi() << endl;
  // cout << "XZ: " << n1.XZGetChisquare() << endl;
  // cout << "YZ: " << n1.YZGetChisquare() << endl;
  // cout << "XY: " << n1.XYGetChisquare() << endl;

  cout << '\n' << "MANUALE" << '\n' << endl;
  cout << "xy: " << n1.YXGetChisquare_m() << endl;
  cout << "xz: " << n1.ZXGetChisquare_m() << endl;
  cout << "yz: " << n1.YZGetChisquare_m() << endl;
  cout << "YXParameter1: " << n1.YXGetParameter(1) << endl;
  cout << "ZXParameter1: " << n1.ZXGetParameter(1) << endl;
  cout << "YZParameter1: " << n1.YZGetParameter(1) << endl;
  cout << "theta: " << n1.GetTheta() << endl;
  cout << "phi: " << n1.GetPhi() << endl;

   // cout << "ALTRA PROVA..." << endl;

   // n1.yzgraph->Fit("1 ++ x");
   // TF1* provafit = n1.yzgraph->GetFunction("1 ++ x");
   // provafit->DrawClone("SAME");

  double par0, par1;
  int valid;
  //USANDO LINEAR LSM
  cout << "Linear LSM..." << endl;
  n1.yzgraph->LeastSquareLinearFit(3,par0,par1,valid,-100,100);
  cout << "par0: " << par0 << endl;
  cout << "par1: " << par1 << endl;
  cout << "valid: " << valid << endl;
//   cout << "USANDO LINEAR FITTER..." << endl;

//   TLinearFitter* yxlf = new TLinearFitter();
//   TLinearFitter* zxlf = new TLinearFitter();
//   TLinearFitter* yzlf = new TLinearFitter();
//   double* x1 = new double;
//   *x1 = 21.7707;
//   double* x2 = new double;
//   *x2 = 31.9766;
//   double* x3 = new double;
//   *x3 = 13.9974;
  
//   TF1* yxprovafit = new TF1("yxprovafit","pol1",-1000,1000);
//   TF1* zxprovafit = new TF1("zxprovafit","pol1",-1000,1000);
//   TF1* yzprovafit = new TF1("yzprovafit","pol1",-1000,1000);

//   yxlf->SetFormula(yxprovafit);
//   zxlf->SetFormula(zxprovafit);
//   yzlf->SetFormula(yzprovafit);

// yzlf->AddPoint(x1, 53.2);
//   yzlf->AddPoint(x2, 0);
//   yzlf->AddPoint(x3, -52.8);

//   yzlf->Eval();

//   cout << "YZGetParameter1: " <<  yzlf->GetParameter(1) << endl;
//   yzlf->DrawClone("same");



  // double chi = n1.XYGetChisquare();
  // cout << "chi square: " << chi << endl;
  //  cout << "Theta: " << n1.GetTheta() << endl;
  // cout << "Phi: " << n1.GetPhi() << endl;
};
  

