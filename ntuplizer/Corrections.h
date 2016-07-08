#include "TF1.h"

static double centMinCoarseC[] = {0,20,60,9999999};

class Corrections{
 public:  
  Corrections(){

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

  float factor(double pt, double eta, int hiBin){
    return f[hiBin]->Eval(pt);
  }
  TF1* f[200];
  TF1* fcent;

};
