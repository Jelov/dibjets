#ifndef CONFIG_H
#define CONFIG_H

#include "TString.h"
#include "parsecode.h"
#include "physics.h"



bool file_exist(const char *fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

class Config {
public:
  TString workdir = "/Users/istaslis/Documents/CMS/bjet2015";

  TString tuplesfolderPbPb = workdir+"/ntuples/eta1p5_jecv2";
  TString tuplesfolderpp = workdir+"/ntuples/eta1p5";

  TString tagcorfolder = workdir+"/dibjets/correctionfiles";


  // void setuptagging(int csvmode)
  // {
  //   if (csvmode==0) {
  //     csvcut = 0.9;
  //     tagcorfolder = workdir+"/dibjets/correctionfiles/tagging_jtsignal2v4_0p90";
  //   }
  //   if (csvmode==1) {
  //     csvcut = 0.85;
  //     tagcorfolder = workdir+"/dibjets/correctionfiles/tagging_jtsignal2v4_0p85";
  //   }
  //   if (csvmode==2) {
  //     csvcut = 0.95;
  //     tagcorfolder = workdir+"/dibjets/correctionfiles/tagging_jtsignal2v4_0p95";
  //   }
  //   cout<<"Using csv = "<<csvcut<<", tagging corrections at "<<tagcorfolder<<endl;
  // }


  TString PbPbjetalgo = "akPu4PF";
  TString ppjetalgo = "ak4PF";
  
  TString getFileName_djt(TString sample)
  {
    TString algo = isPbPb(sample) ? PbPbjetalgo : ppjetalgo;
    TString folder = isPbPb(sample) ? tuplesfolderPbPb : tuplesfolderpp;
    TString fname = folder+   // ((mc(sample)  && isPbPb(sample)) ? "/privmc" : "") 
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
    TString folder = isPbPb(sample) ? tuplesfolderPbPb : tuplesfolderpp;    
    TString fname = folder+"/"+sample+algo+"_inc.root";
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