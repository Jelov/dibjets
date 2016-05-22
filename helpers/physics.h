#ifndef PHYSICS_H
#define PHYSICS_H

#include "TString.h"
#include "parsecode.h"


//PHYSICAL CONSTANTS

const float pthatcut = 50;

const float pt1cut = 100;
const float pt2cut = 40;

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
vector<float> processWeights = {1.2,1.,0.04,0.04};
// vector<float> processWeights = {1.232,1.,0.266,0.04};

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
vector<float> bkgfractionInNearSide = {0.8726,0.4479,0.0801};

//true
// vector<float> bkgfractionInNearSide = {0.8455,0.3612,0.0311};

// vector<float> bkgfractionInNearSide = {0.96,0.55,0.18};



bool applyTaggingCorrection = true;
bool applyTriggerCorr = true;



#endif