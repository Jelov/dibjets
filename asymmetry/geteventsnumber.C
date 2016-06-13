#include "../helpers/config.h"

TTree *nt;

vector<TString> trigs;

void getentries(TString cut, bool showxj = false)
{
  cout<<cut<<endl;
  auto wh = new TH1F("wxj","wxj",100,0,1);

  for (auto t:trigs) {
    cout<<t<<endl;

    if (t=="all events") t="1";

     //not correct: events with 0 weight can be repeated events, e.g. triggerPt>80 but CaloJet60
     // nt->Project("xj","jtpt2/jtpt1",Form("%s && %s",t.Data(),cut.Data()));
     //float r= h->Integral();  

    float r = nt->Project("wxj","jtpt2/jtpt1",Form("weight*(%s && %s)",t.Data(),cut.Data()));
    float e = wh->Integral();

    cout<<"   raw: "<<r<<" weighted: "<<e;
    if (showxj) 
      cout<<" xJ = "<<wh->GetMean()<<"±"<<wh->GetMeanError()<<endl;
    else cout<<endl;
  }

  cout<<endl;
  delete wh;
}

void getnumbers(TString file)
{
  auto f = config.getfile_djt(file);
  nt = (TTree *)f->Get("nt");
  getentries("bin<20 && numTagged<=6");
  getentries("bin<20 && numTagged<=6 && jtpt1>100");
  getentries("bin<20 && numTagged<=6 && jtpt1>100 && jtpt2>40 && dphi21>2.1");
  getentries("bin<20 && numTagged<=6 && jtpt1>100 && jtpt2>40 && dphi21>2.1 && discr_csvV1_1>0.9 && discr_csvV1_2>0.9",true);

}

void geteventsnumber()
{
  trigs = {"all events","hltCSV80","hltCSV60 && !hltCSV80"};
  cout<<"-------CSV triggers-------"<<endl;
  getnumbers("dtPbbjt");

  trigs = {"all events","hltCaloJet100","hltCaloJet80 && !hltCaloJet100","hltCaloJet60 && !hltCaloJet80 && !hltCaloJet100"};
  cout<<"-------Calo jet triggers-------"<<endl;
  getnumbers("dtPbjcl");
}