#ifndef PHYSICS_H
#define PHYSICS_H

#include "TString.h"
#include "parsecode.h"
#include "looptuple.h"


//PHYSICAL CONSTANTS

const float pthatcut = 50;

const float pt1cut = 100;
const float pt2cut = 40;

const float etacut = 1.5;

float csvcut1 = 0.9; //modifiable later
float csvcut2 = 0.9; //modifiable later

const float PI = 3.141593;
const float PI23 = PI*2/3;
const float PI13 = PI*1/3;

const float NaN = -999;


//DEFINITION OF TAGGER AND PARTNER JET
TString discr_csvV1_1 = "discr_csvV1_1";
TString discr_csvV1_2 = "discr_csvV1_2";
TString discr_csvV1_Signal2 = "discr_csvV1_Signal2";
TString jtptSL = "jtptSL";
TString dphiSL1 = "dphiSL1";
TString jtetaSL = "jtetaSL";
TString subidSL = "subidSL";
TString refptSL = "refptSL";
TString pairCodeSL1 = "pairCodeSL1";
TString refparton_flavorForBSL = "refparton_flavorForBSL";
TString SLord = "SLord";



//process code reweighting: GSP, FCR, FEX, FEX2
// vector<float> processWeights = {1.2,1.,0.04,0.04};
// vector<float> processWeights = {1.2,1.,0.04,1.};
vector<float> processWeights = {1.089,1.,0.216,1};


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


//definition of the embedded signal
bool IsSignal(dict d) { return d[subidSL]==0 && d[refptSL]>20;}

bool IsSignal(float subid, float refpt) {return subid==0 && refpt>20;}


//physical algo shortcuts
float weight1SLpp(dict d)
{
  float w = d["weight"];
  if (d[pairCodeSL1]==0) w*=processweight((int)d["bProdCode"]);
  return w;
}

float weight1SLPbPb(dict d)
{
  float w = d["weight"];
  if (d[pairCodeSL1]==0 && IsSignal(d)) w*=processweight((int)d["bProdCode"]);
  return w;
}

float weight12(dict d)
{
  float w = d["weight"];
  if (d["pairCode21"]==0 && IsSignal(d["subid2"],d["refpt2"])) w*=processweight((int)d["bProdCode"]);
  return w;
}

float weight1Signal2(dict d)
{
  float w = d["weight"];
  if (d["pairCodeSignal21"]==0 && IsSignal(d["subidSignal2"],d["refptSignal2"])) w*=processweight((int)d["bProdCode"]);
  return w;
}

float weight(float weight, float pairCode, float subid, float refpt, float bProdCode)
{
  float w = weight;
  if ((int)pairCode==0 && IsSignal(subid,refpt)) w*=processweight((int)bProdCode);
  return w;
}


bool NearSide(dict d)
{
  //return d["dphiSL1"]<PI13;

  float dphi = d[dphiSL1];
  float deta = abs(d["jteta1"]-d[jtetaSL]); 
  return (dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1);
}

bool AwaySide(dict d)
{
  return d[dphiSL1]>PI23;
}


//for pthat>50, fullMC
// vector<float> bkgfractionInNearSide = {0.8692479730,0.4500232041,0.0801286325};
//|eta|<2
// vector<float> bkgfractionInNearSide = {0.8726,0.4479,0.0801};

//|eta|<1.5
vector<float> bkgfractionInNearSide = {0.7548671365,0.3059994578,0.0966218933};


bool applyTaggingCorrection = true;
bool applyTriggerCorr = true;



#endif