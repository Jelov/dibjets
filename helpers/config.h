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

bool file_exist(const char *fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

class Config {
public:
  TString tuplesfolder = "/Users/istaslis/Documents/CMS/bjet2015/ntuples";//"/data_CMS/cms/lisniak/bjet2015";
  TString PbPbjetalgo = "akPu4PF";
  TString ppjetalgo = "ak4PF";
  
  TString getFileName_djt(TString sample)
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString fname = tuplesfolder+    //(dt(sample) ? "/latestdata" : "") 
                      +"/"+sample+algo+"_djt.root";
    if (!file_exist(fname))
      cout<<"File "<<fname<<" doesn\'t exist"<<endl;
    return fname;
  }

  TFile *getfile_djt(TString sample)
  {
    return new TFile(getFileName_djt(sample));
  }

  TString getFileName_inc(TString sample)
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString fname = tuplesfolder+"/"+sample+algo+"_inc.root";
    if (!file_exist(fname))
      cout<<"File "<<fname<<" doesn\'t exist"<<endl;
    return fname;
  }
};

Config config;