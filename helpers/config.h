#include "TString.h"
#include "parsecode.h"

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
  
  TString getFileName_djt(TString sample, TString cbin = "")
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString cbinf = cbin!="" ? "cbin"+cbin+"/" : "";
    TString fname = tuplesfolder+    //(dt(sample) ? "/latestdata" : "") 
                      +"/"+cbinf+sample+algo+"_djt.root";
    if (!file_exist(fname))
      cout<<"File "<<fname<<" doesn\'t exist"<<endl;
    return fname;
  }
  TString getFileName_inc(TString sample, TString cbin = "")
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString cbinf = cbin!="" ? "cbin"+cbin+"/" : "";
    TString fname = tuplesfolder+"/"+cbinf+sample+algo+"_inc.root";
    if (!file_exist(fname))
      cout<<"File "<<fname<<" doesn\'t exist"<<endl;
    return fname;
  }
};

Config config;