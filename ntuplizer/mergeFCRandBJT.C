#include <fstream>
#include "../helpers/parsecode.h"
#include "TChain.h"

int getpthatind(vector<int> &vpthat, float pthat)
{
  for (unsigned i=0;i<vpthat.size();i++)
    if (pthat<vpthat[i]) return i-1;
  return vpthat.size()-1;
}

void mergeFCRandBJT(TString fcrsample, TString bjtsample, TString outsample)
{
  cout<<"Merging b-jet FCR and filtered samples..."<<endl;

  vector<float> fcrCS = {1.890e-03,1.242e-04,1.348e-05,4.015E-06,1.468e-06,4.831e-07,1.889e-07};
  vector<int> pthats = { 15,   30,   50,   65,   80,   100,     120};


  TString djtfcrsample(fcrsample);
  djtfcrsample.ReplaceAll("_inc","_djt");
  TString djtbjtsample(bjtsample);
  djtbjtsample.ReplaceAll("_inc","_djt");

  //calcualte Nev with FCR
  vector<float> w (pthats.size());
  auto ffcr = new TFile(djtfcrsample);
  auto ntfcr = (TTree *)ffcr->Get("nt");

  int nFCR = ntfcr->GetEntries();

  auto fbjt = new TFile(djtbjtsample);
  auto ntbjt = (TTree *)fbjt->Get("nt"); 

  for (unsigned i=0;i<pthats.size();i++) {
    int Nfcr = 0 , Nbjt = 0, Nbjtfc = 0;

    int pt1 = pthats[i];
    int pt2 = i<pthats.size()-1 ? pthats[i+1] : 1E6;

    float cs = i<pthats.size()-1 ? fcrCS[i]-fcrCS[i+1]:fcrCS[i];


    Nfcr=ntfcr->GetEntries(Form("pthat>%d && pthat<=%d",pt1,pt2));
    Nbjtfc=ntbjt->GetEntries(Form("pthat>%d && pthat<=%d && bProdCode==1",pt1,pt2));
    cout<<"Nfcr "<<Nfcr<<" Nbjtfc "<<Nbjtfc<<endl;
    w[i] = cs/(float)(Nfcr+Nbjtfc);
    //w[i] = Nbjtfc/(float)(Nfcr+Nbjtfc);

    //cross-check:
    //Nbjt=ntbjt->GetEntries(Form("pthat>%d && pthat<=%d",pt1,pt2));
    //if (i<pthats.size()-1) 
    //  cout<<fcrCS[i]-fcrCS[i+1]<<" = "<<(CS[i]-CS[i+1])*filterefficiency[i]*Nbjtfc/Nbjt<<endl;


  }

  delete ntfcr;
  delete ntbjt;

  cout<<"Weights for fcr ";for (auto ww:w) cout<<ww<<" "; cout<<endl;


  TChain ch("nt");
  ch.Add(fcrsample);
  ch.Add(bjtsample);
  ch.Merge(outsample);

  //now I loop and update if it's fcr


  auto f = new TFile(outsample,"update");
  auto nt = (TTree *)f->Get("nt");

  float fcrweight, pthatweight, bProdCode, pthatsample,pthat;

  TBranch *bw;

  bw =  nt->Branch("fcrweight",&fcrweight);
  nt->SetBranchAddress("pthatweight",&pthatweight);
  nt->SetBranchAddress("bProdCode",&bProdCode);
  nt->SetBranchAddress("pthatsample",&pthatsample);
  nt->SetBranchAddress("pthat",&pthat);

  int n = nt->GetEntries();

  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;

    nt->GetEntry(i);

    fcrweight = 1;
    if (i<nFCR || bProdCode==1) {
      int ind = getpthatind(pthats,pthat);
      if (ind!=-1) fcrweight = w[ind];
    } else fcrweight = pthatweight;

    if (fcrweight<=0) 
      cout<<pthatweight<<" "<<fcrweight<<" "<<pthatsample<<" "<<i<<endl;


    bw->Fill();
  }


  nt->Write("nt",TObject::kOverwrite);

  f->Close();


}

void mergeFCRandBJT()
{
  TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
  TString jetalgo = "ak4PF";

  // mergeFCRandBJT(outputfolder+"mcppbfc"+jetalgo+"_djt.root",
		//  outputfolder+"mcppbjt"+jetalgo+"_djt.root",
		//  outputfolder+"mcppbfa"+jetalgo+"_djt.root");


  mergeFCRandBJT(outputfolder+"/mcPbbfcakVs4PF_inc.root", outputfolder+"/mcPbbjtakVs4PF_inc.root", outputfolder+"/mcPbbfaakVs4PF_inc.root");
 
}
