#include "TF1.h"

static double centMinCoarseC[] = {0,20,60,9999999};

class Parametrization{
 public:  
  Parametrization(){

    fcent = new TF1("fcent","expo",0,200);
    fcent->SetParameters(5.63,-0.0455);

    for(int i = 0; i < 200; ++i){
      f[i] = new TF1(Form("f%d",i),"1./([0]-[1]*exp(-[2]*x)+[3]/(x-[4])/(x-[4]))",0,1000);
      f[i]->SetParameters(1.,0.08,0.01,fcent->Eval(i),15);
    }    
  }
  
  int coarseCentrality(int hiBin){
    int c = 0;
    while(hiBin > centMinCoarseC[c+1]) c++;
    return c;
  }

  double factor(double pt, double eta, int hiBin){
    return f[hiBin]->Eval(pt);
  }
  TF1* f[200];
  TF1* fcent;

};

class Corrections : public Parametrization{};


class Smearing : public Parametrization{
 public:
  Smearing(double c=0, double s=0, double nmax = 16, double nslope = -0.1){

    gaus = new TF1("fgaus","gaus",-10,10);
    gaus->SetParameters(1,0,1);
    fcent = new TF1("freso","[0]+[1]*x",0,200);
    fcent->SetParameters(nmax,nslope);

    for(int i = 0; i < 200; ++i){
      f[i] = new TF1(Form("freso%d",i),"sqrt(pow([0],2)+pow([1]/sqrt(x),2)+pow([2]/x,2))",0,1000);

      double r = fcent->Eval(i);
      if(r < 0) r = 0;
      f[i]->SetParameters(c,s,fcent->Eval(i));
    }


  }

  double roll(double pt, double eta, int hiBin){
    return gaus->GetRandom()*factor(pt,eta,hiBin);
  }

  double rollpp(int hiBin){
    return gaus->GetRandom()*fcent->Eval(hiBin);
  }


  TF1* gaus;

};
