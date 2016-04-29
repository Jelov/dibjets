#ifndef PHYSICS_H
#define PHYSICS_H

#include "TString.h"
#include "parsecode.h"


//PHYSICAL CONSTANTS
const float pt1cut = 100;
const float pt2cut = 40;

const float PI = 3.141593;
const float PI23 = PI*2/3;
const float PI13 = PI*1/3;

const float NaN = -999;

//process code reweighting: GSP, FCR, FEX, FEX2
vector<float> processWeights = {1.2,1.,0.04,0.04};
float processweight(int bProdCode)
{
  if (bProdCode>=0 && bProdCode<=4)
    return processWeights[bProdCode];
  else {
    cout<<"Process code "<<bProdCode<<" is wrong!!!!!"<<endl;
    return NaN;
  }
}

//centrality bins
int Nbins = 3;
vector<float> bins = {0,20,60,200};
vector<TString> binnames = {"0-10%", "10-30%", "30-100%"}; 

int getbinindex(float bin)
{
  for(unsigned i=0;i<bins.size();i++)
    if (bins[i]>bin) return i-1;
  return bins.size()-2;
}


//physical algo shortcuts
float weight1SL(dict d)
{
  float w = d["weight"];
  if (d["pairCodeSL1"]==0) w*=processweight((int)d["bProdCode"]);
  return w;
}

//definition of the embedded signal
bool IsSignal(dict d) { return d["subidSL"]==0 && d["refptSL"]>20;}

bool NearSide(dict d)
{
  //return d["dphiSL1"]<PI13;

  float dphi = d["dphiSL1"];
  float deta = abs(d["jteta1"]-d["jtetaSL"]); 
  return (dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1);
}

bool AwaySide(dict d)
{
  return d["dphiSL1"]>PI23;
}


//Hydjet/(Hydjet+Signal) coefficients for NearSide
//check out mistag/hydjetestimation.C and mistag/hydjetclosure.C for more information
//fancy subtraction
vector<float> bkgfractionInNearSide = {0.8690670729,0.4207638502,0.0000000000};
//simple dphi
//vector<float> bkgfractionInNearSide = {0.8502221704,0.4357484281,0.010240517};



#endif