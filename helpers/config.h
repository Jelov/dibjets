#ifndef CONFIG_H
#define CONFIG_H

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
  
  TString getFileName_djt(TString sample)
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString fname = tuplesfolder+   // ((mc(sample)  && isPbPb(sample)) ? "/privmc" : "") 
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

  TFile *getfile_inc(TString sample)
  {
    return new TFile(getFileName_inc(sample));
  }
};

Config config;

#endif