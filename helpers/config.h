#ifndef CONFIG_H
#define CONFIG_H

#include "./CMS_lumi.C"
#include "TString.h"
#include "parsecode.h"
// #include "physics.h"

int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

// second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
int iPos=11; // : top-left, left-aligned
//int  iPos=33;// : top-right, right-aligned
// iPos=22 : center, centered
// mode generally : 
//   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

TString lumi_sqrtSpp = "25.8 pb^{-1} (5.02 TeV pp)";
TString lumi_sqrtSPbPb = "404 #mub^{-1} (5.02 TeV PbPb)";


bool file_exist(const char *fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

bool firstRunMacro = true;

class Config {
public:
  TString workdir = "/Users/istaslis/Documents/CMS/bjet2015";

  TString tuplesfolderPbPb = workdir+"/ntuples/eta1p5_jecv2";
  TString tuplesfolderpp = workdir+"/ntuples/eta1p5";

  // TString tuplesfolderPbPb = workdir+"/ntuples/eta1p5";
  // TString tuplesfolderpp = workdir+"/ntuples/eta1p5";

  TString tagcorfolder = workdir+"/dibjets/correctionfiles";

  TString PbPbjetalgo = "akPu4PF";
  TString ppjetalgo = "ak4PF";
  
  Config()
  {

    writeExtraText = true;       // if extra text
    extraText  = "Preliminary";  // default extra text is "Preliminary"
    //lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    lumi_sqrtS = lumi_sqrtSPbPb;       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)


  }


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