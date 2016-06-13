#include <iostream>
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "../helpers/config.h"

using namespace std;


Double_t fit_pt(Double_t *x, Double_t *par)
{
  Double_t xv =x[0];
  Double_t f = (xv>=100. && xv<=160.)*(par[0]+xv*par[1]+xv*xv*par[2]) + (xv>160.&&xv<=200.)*( par[0]+160.*par[1]+160.*160.*par[2]-160.*par[3]-160.*160.*par[4] + xv*par[3] + xv*xv*par[4]);
  return f;
}

Double_t fit_pt2(Double_t *x, Double_t *par)
{
  Double_t xv =x[0];
  Double_t f = (xv>=100. && xv<=160.)*(par[0]+xv*par[1]) + (xv>160.&&xv<=200.)*( par[0]+160.*par[1]-160.*par[2]-160.*160.*par[3] + xv*par[2] + xv*xv*par[3]);
  return f;
}

Double_t fit_pt3(Double_t *x, Double_t *par)
{
  Double_t xv =x[0];
  Double_t f = (xv>=100. && xv<=160.)*(par[0]+xv*par[1]) + (xv>160.&&xv<=200.)*( par[0]+160.*par[1]-160.*par[2]-160.*160.*par[3] - 160.*160.*160.*par[4]+ xv*par[2] + xv*xv*par[3] + xv*xv*xv*par[4]);
  return f;
}


void deriveOffTagEff(bool ppPbPb=false, bool requireFCR=false, bool doRebin= true, bool useGenPt=false, int cBin=1, bool doResid=false, bool savePlots=true)
{

  // cBin =0 --> 1-100%
  // cBin =1 --> 1-10%
  // cBin =2 --> 10-30%
  // cBin =3 --> 30-100%

  if(!ppPbPb) cBin=0;
  
   TH1::SetDefaultSumw2();

  TFile *fdjt = NULL;

  if(ppPbPb)  fdjt = config.getfile_djt("mcPbbfa");
  else fdjt = config.getfile_djt("mcppbfa");

  TNtuple *nt = (TNtuple *)fdjt->Get("nt");

  float bProdCode, bin, refparton_flavorForB1,refparton_flavorForBSB;
  float weight, pthat, jtpt1, refpt1, jteta1, discr_csvV1_1, jtptSB, refptSB, jtetaSB, discr_csvV1_SB, dphiSB1;
  
  nt->SetBranchAddress("bProdCode",&bProdCode);
  nt->SetBranchAddress("bin",&bin);
   nt->SetBranchAddress("weight",&weight);
   nt->SetBranchAddress("pthat",&pthat);
  nt->SetBranchAddress("jtpt1",&jtpt1);
  nt->SetBranchAddress("refpt1",&refpt1);
  nt->SetBranchAddress("jteta1",&jteta1);
  nt->SetBranchAddress("discr_csvV1_1",&discr_csvV1_1);
  nt->SetBranchAddress("refparton_flavorForB1",&refparton_flavorForB1);
  nt->SetBranchAddress("jtptSB",&jtptSB);
  nt->SetBranchAddress("refptSB",&refptSB);
  nt->SetBranchAddress("jtetaSB",&jtetaSB);
  nt->SetBranchAddress("discr_csvV1_SB",&discr_csvV1_SB);
  nt->SetBranchAddress("refparton_flavorForBSB",&refparton_flavorForBSB);
  nt->SetBranchAddress("dphiSB1",&dphiSB1);
  

  float maxPtLead = 200;
  float maxPtPart = 140;
  
  TH1F *hLeadPt=new TH1F("hLeadPt","hLeadPt;Leading jet p_{T} (GeV);Counts",20,100,maxPtLead);
  TH1F *hLeadEta=new TH1F("hLeadEta","hLeadEta;Leading jet #eta;Counts",20,-2.,2.);
  TH1F *hPartPt=new TH1F("hPartPt","hPartPt;Partner jet p_{T} (GeV);Counts",20,40,maxPtPart);
  TH1F *hPartEta=new TH1F("hPartEta","hPartEta;Partner jet #eta;Counts",20,-2.,2.);
  TH1F *hCent=new TH1F("hCent","hCent;Centrality (%);Counts",10,0,100);
  TH1F *hXj=new TH1F("hXj","hXj;x_{J};Counts",10,0,1);

  
  TH1F *hTagLeadPt=new TH1F("hTagLeadPt","hTagLeadPt;Leading jet p_{T} (GeV);Counts",20,100,maxPtLead);
  TH1F *hTagLeadEta=new TH1F("hTagLeadEta","hTagLeadEta;Leading jet #eta;Counts",20,-2.,2.);
  TH1F *hTagPartPt=new TH1F("hTagPartPt","hTagPartPt;Partner jet p_{T} (GeV);Counts",20,40,maxPtPart);
  TH1F *hTagPartEta=new TH1F("hTagPartEta","hTagPartEta;Partner jet #eta;Counts",20,-2.,2.);
  TH1F *hTagLeadCent=new TH1F("hTagLeadCent","hTagLeadCent;Centrality (%);Counts",10,0,100);
  TH1F *hTagPartCent=new TH1F("hTagPartCent","hTagPartCent;Centrality (%);Counts",10,0,100);
  
  TH1F *hTagEffLeadPt=new TH1F("hTagEffLeadPt","hTagEffLeadPt;Leading jet p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hTagEffLeadEta=new TH1F("hTagEffLeadEta","hTagEffLeadEta;Leading jet #eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hTagEffPartPt=new TH1F("hTagEffPartPt","hTagEffPartPt;Partner jet p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hTagEffPartEta=new TH1F("hTagEffPartEta","hTagEffPartEta;Partner jet #eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hTagEffLeadCent=new TH1F("hTagEffLeadCent","hTagEffLeadCent;Centrality (%);Tagging Efficiency",10,0,100);
  TH1F *hTagEffPartCent=new TH1F("hTagEffPartCent","hTagEffPartCent;Centrality (%);Tagging Efficiency",10,0,100);
  TH1F *hTagXj=new TH1F("hTagXj","hXj;x_{J};Counts",10,0,1);

  
  TH1F *hw1TagLeadPt=new TH1F("hw1TagLeadPt","hw1TagLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw1TagLeadEta=new TH1F("hw1TagLeadEta","hw1TagLeadEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw1TagPartPt=new TH1F("hw1TagPartPt","hw1TagPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw1TagPartEta=new TH1F("hw1TagPartEta","hw1TagPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw1TagCent=new TH1F("hw1TagCent","hw1TagCent;p_{T} (GeV);Tagging Efficiency",10,0,100);
  TH1F *hw1TagXj=new TH1F("hw1TagXj","hXj;x_{J};Counts",10,0,1);
  
  TH1F *hw1TagEffLeadPt=new TH1F("hw1TagEffLeadPt","hw1TagEffLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw1TagEffLeadEta=new TH1F("hw1TagEffLeadEta","hw1TagEffLeadEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw1TagEffPartPt=new TH1F("hw1TagEffPartPt","hw1TagEffPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw1TagEffPartEta=new TH1F("hw1TagEffPartEta","hw1TagEffPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw1TagEffPartCent=new TH1F("hw1TagEffPartCent","hw1TagEffPartCent;p_{T} (GeV);Tagging Efficiency",10,0,100);
  TH1F *hw1TagEffXj=new TH1F("hw1TagEffXj","hXj;x_{J};Counts",10,0,1);
  
  TH1F *hw2TagLeadPt=new TH1F("hw2TagLeadPt","hw2TagLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw2TagLeadEta=new TH1F("hw2TagLeadEta","hw2TagLeadEta;#eta;Tagging Efficiency",20,-2.,2.);  
  TH1F *hw2TagPartPt=new TH1F("hw2TagPartPt","hw2TagPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw2TagPartEta=new TH1F("hw2TagPartEta","hw2TagPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw2TagCent=new TH1F("hw2TagCent","hw2TagCent;p_{T} (GeV);Tagging Efficiency",10,0,100);

  TH1F *hw2TagEffLeadPt=new TH1F("hw2TagEffLeadPt","hw2TagEffLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw2TagEffLeadEta=new TH1F("hw2TagEffLeadEta","hw2TagEffLeadEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw2TagEffPartPt=new TH1F("hw2TagEffPartPt","hw2TagEffPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw2TagEffPartEta=new TH1F("hw2TagEffPartEta","hw2TagEffPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw2TagEffCent=new TH1F("hw2TagEffCent","hw2TagEffCent;p_{T} (GeV);Tagging Efficiency",10,0,100);
  
  TH1F *hw3TagLeadPt=new TH1F("hw3TagLeadPt","hw3TagLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw3TagLeadEta=new TH1F("hw3TagLeadEta","hw3TagLeadEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw3TagPartPt=new TH1F("hw3TagPartPt","hw3TagPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw3TagPartEta=new TH1F("hw3TagPartEta","hw3TagPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw3TagCent=new TH1F("hw3TagCent","hw3TagCent;p_{T} (GeV);Tagging Efficiency",10,0,100);

  TH1F *hw3TagEffLeadPt=new TH1F("hw3TagEffLeadPt","hw3TagEffLeadPt;p_{T} (GeV);Tagging Efficiency",20,100,maxPtLead);
  TH1F *hw3TagEffLeadEta=new TH1F("hw3TagEffLeadEta","hw3TagEffLeadEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw3TagEffPartPt=new TH1F("hw3TagEffPartPt","hw3TagEffPartPt;p_{T} (GeV);Tagging Efficiency",20,40,maxPtPart);
  TH1F *hw3TagEffPartEta=new TH1F("hw3TagEffPartEta","hw3TagEffPartEta;#eta;Tagging Efficiency",20,-2.,2.);
  TH1F *hw3TagEffCent=new TH1F("hw3TagEffCent","hw3TagEffCent;p_{T} (GeV);Tagging Efficiency",10,0,100);
  
  
  
  Long64_t nentries = nt->GetEntries();

  double totCounts0 =0;

  // loop 0:  Fill centrality.  PbPb only

  TF1 *fLeadCent = new TF1("fLeadCent","pol3",0,100);
  TF1 *fPartCent = new TF1("fPartCent","pol3",0,100);
  //  TF1 *fLeadPt = new TF1("fLeadPt","pol4",100,maxPtLead);
  TF1 *fLeadPt = new TF1("piecewisePt1",fit_pt2,100,maxPtLead,4);
  //  TF1 *fLeadPt = new TF1("piecewisePt1",fit_pt3,100,maxPtLead,5);
  TF1 *fLeadEta = NULL;
  //fLeadEta=new TF1("fLeadEta","pol4",-2.,2.);
  fLeadEta=new TF1("fLeadEta","[0]+[1]*x*x+[2]*x*x*x*x",-2.,2.);
  //else fLeadEta=new TF1("fLeadEta","pol6",-2.,2.);
  TF1 *fPartEta = NULL;
  //fPartEta=new TF1("fPartEta","pol4",-2.,2.);
  fPartEta=new TF1("fPartEta","[0]+[1]*x*x+[2]*x*x*x*x",-2.,2.);
  //else fPartEta=new TF1("fPartEta","pol6",-2.,2.);
  //TF1 *fPartPt = new TF1("fPartPt","pol4",40,maxPtPart);
  TF1 *fPartPt = NULL;
  //if(ppPbPb)fPartPt=new TF1("fPartPt","pol5",40,maxPtPart);
  //else fPartPt=new TF1("fPartPt","pol3",40,maxPtPart);
  if(ppPbPb&&cBin==3)fPartPt=new TF1("fPartPt","pol2",40,maxPtPart);
  else fPartPt=new TF1("fPartPt","pol3",40,maxPtPart);
  TF1 *fLeadPt2 = new TF1("fLeadPt2","gaus",100,maxPtLead);
  TF1 *fPartPt2=new TF1("fPartPt2","pol5",40,maxPtPart);



  // loop 1:  Leading jet efficiency
  // centrality
  
  if(ppPbPb){

    for (Long64_t i=0; i<nentries;i++) {
      nt->GetEntry(i);
      
      if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
      else if(!requireFCR && bProdCode==2) continue;
      
      float pt1=jtpt1;
      if(useGenPt) pt1 = refpt1;
      
      if(pt1<100) continue;
      if(abs(refparton_flavorForB1)!=5) continue;
      
      float ptSB=jtptSB;
      if(useGenPt) ptSB = refptSB;
      
      if(ptSB<40) continue;
      if(dphiSB1<2./3.*acos(-1.)) continue;
      if(refptSB<20) continue;
	   
      hCent->Fill(bin/2.,weight);
      
      totCounts0 += weight;
      
      if(discr_csvV1_1<0.9) continue;
      //if(discr_csvV1_SB<0.9) continue;
      
      hTagLeadCent->Fill(bin/2.,weight);
      
    }

    
    hTagEffLeadCent->Divide(hTagLeadCent,hCent,1.,1.,"B");  
    hTagEffLeadCent->Fit(fLeadCent);
  }

  cout<<" totCounts0 "<<totCounts0<<endl;
  
  // loop 2:  Leading jet efficiency
  // pt & eta

  
  
  double totCounts0a=0;
  double totCounts1=0;
  double totCorr1=0;
  
  for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
    
     if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
     else if(!requireFCR && bProdCode==2) continue;

       if(cBin==1){
	 if(bin>=20) continue;
       }
       else if(cBin==2){
	 if(bin<20 || bin>=60) continue;       
     }
       else if(cBin==3){
	 if(bin<60) continue;       
     }
     
     float pt1=jtpt1;
     if(useGenPt) pt1 = refpt1;
     
     if(pt1<100) continue;
     if(abs(refparton_flavorForB1)!=5) continue;

     float ptSB=jtptSB;
     if(useGenPt) ptSB = refptSB;
     
     if(ptSB<40) continue;
     if(dphiSB1<2./3.*acos(-1.)) continue;
     if(refptSB<20) continue;
      
     hLeadPt->Fill(pt1,weight);
     hLeadEta->Fill(jteta1,weight);
     hPartPt->Fill(ptSB,weight);
     hPartEta->Fill(jtetaSB,weight);
     hXj->Fill(ptSB/pt1,weight);
     
     totCounts0a += weight;
     
     if(discr_csvV1_1<0.9) continue;
	      
     float centTagEffCorr = 1.;
     if(ppPbPb) centTagEffCorr=1./fLeadCent->Eval(bin/2.);

     totCounts1 += weight;
     totCorr1+=weight*centTagEffCorr;
     
     hTagLeadPt->Fill(pt1,weight*centTagEffCorr);
     hTagLeadEta->Fill(jteta1,weight*centTagEffCorr);
     
     
  }

  cout<<" totCounts0a "<<totCounts0a<<endl;
  cout<<" totCounts1 "<<totCounts1<<endl;
  
  hTagLeadPt->Scale(totCounts1/totCorr1);
  hTagLeadEta->Scale(totCounts1/totCorr1);
 
  hTagEffLeadPt->Divide(hTagLeadPt,hLeadPt,1.,1.,"B");  
  hTagEffLeadEta->Divide(hTagLeadEta,hLeadEta,1.,1.,"B");  


  hTagEffLeadPt->Fit(fLeadPt);
  hTagEffLeadEta->Fit(fLeadEta);


 // loop 3:  Partner jet efficiency
  // centrality
  
  if(ppPbPb){

    for (Long64_t i=0; i<nentries;i++) {
      nt->GetEntry(i);
      
      if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
      else if(!requireFCR && bProdCode==2) continue;
      
      float pt1=jtpt1;
      if(useGenPt) pt1 = refpt1;
      
      if(pt1<100) continue;
      if(abs(refparton_flavorForB1)!=5) continue;
      
      float ptSB=jtptSB;
      if(useGenPt) ptSB = refptSB;
      
      if(ptSB<40) continue;
      if(dphiSB1<2./3.*acos(-1.)) continue;
      if(refptSB<20) continue;
      
      //hCent->Fill(bin/2.,weight);
      
      totCounts0 += weight;
      
      //if(discr_csvV1_1<0.9) continue;
      if(discr_csvV1_SB<0.9) continue;
      
      hTagPartCent->Fill(bin/2.,weight);
      
    }

    
    hTagEffPartCent->Divide(hTagPartCent,hCent,1.,1.,"B");  
    hTagEffPartCent->Fit(fPartCent);
  }
 
 // loop 4:  Partner jet efficiency
  // pt & eta

  
  double totCounts2 =0;
  double totCorr2 =0;
  
 for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
    
     if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
     else if(!requireFCR && bProdCode==2) continue;

       if(cBin==1){
	 if(bin>=20) continue;
       }
       else if(cBin==2){
	 if(bin<20 || bin>=60) continue;       
     }
       else if(cBin==3){
	 if(bin<60) continue;       
     }

     
     float pt1=jtpt1;
     if(useGenPt) pt1 = refpt1;
     
     if(pt1<100) continue;
     if(abs(refparton_flavorForB1)!=5) continue;

     float ptSB=jtptSB;
     if(useGenPt) ptSB = refptSB;
     
     if(ptSB<40) continue;
     if(dphiSB1<2./3.*acos(-1.)) continue;
     if(refptSB<20) continue;
     
     if(discr_csvV1_SB<0.9) continue;


     float centTagEffCorr = 1.;
     if(ppPbPb) centTagEffCorr=1./fPartCent->Eval(bin/2.);

     totCounts2 += weight;
     totCorr2+=weight*centTagEffCorr;
     
     hTagPartPt->Fill(ptSB,weight*centTagEffCorr);
     hTagPartEta->Fill(jtetaSB,weight*centTagEffCorr);

     
     totCounts2+=weight;     
     totCorr2+=weight*centTagEffCorr;

  }

 cout<<" totCorr2 / totCounts2 "<<totCorr2/totCounts2<<endl;

  hTagPartPt->Scale(totCounts2/totCorr2);
  hTagPartEta->Scale(totCounts2/totCorr2);
 
  hTagEffPartPt->Divide(hTagPartPt,hPartPt,1.,1.,"B");  
  hTagEffPartEta->Divide(hTagPartEta,hPartEta,1.,1.,"B");  


  hTagEffPartPt->Fit(fPartPt);
  hTagEffPartEta->Fit(fPartEta);

 // loop 5: Closure test

  double totCounts3 =0;
  double totCorr3 =0;

  
for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
    
     if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
     else if(!requireFCR && bProdCode==2) continue;

       if(cBin==1){
	 if(bin>=20) continue;
       }
       else if(cBin==2){
	 if(bin<20 || bin>=60) continue;       
     }
       else if(cBin==3){
	 if(bin<60) continue;       
     }

     
     float pt1=jtpt1;
     if(useGenPt) pt1 = refpt1;

     if(pt1<100) continue;
     if(abs(refparton_flavorForB1)!=5) continue;
     
     float ptSB=jtptSB;
     if(useGenPt) ptSB = refptSB;

     if(ptSB<40) continue;
     if(dphiSB1<2./3.*acos(-1.)) continue;
     if(refptSB<20) continue;
     
     if(discr_csvV1_1<0.9) continue;
     if(discr_csvV1_SB<0.9) continue;

     float centTagEffCorr = 1.;
     if(ppPbPb) centTagEffCorr=1./fLeadCent->Eval(bin/2.)/fPartCent->Eval(bin/2.);
     
     float leadTagEffPt = fLeadPt->Eval(pt1);
     if(pt1>maxPtLead) leadTagEffPt=fLeadPt->Eval(maxPtLead);
     float leadTagEffEta = fLeadEta->Eval(jteta1);
     float leadTagEffCorr = 1./leadTagEffPt/leadTagEffEta;

     float partTagEffPt = fPartPt->Eval(ptSB);
     if(ptSB>maxPtPart) partTagEffPt=fPartPt->Eval(maxPtPart);
     float partTagEffEta = fPartEta->Eval(jtetaSB);
     float partTagEffCorr = 1./partTagEffPt/partTagEffEta;

     float combTagEffCorr = centTagEffCorr*leadTagEffCorr*partTagEffCorr;
     

     hw1TagLeadPt->Fill(pt1,weight*combTagEffCorr);
     hw1TagLeadEta->Fill(jteta1,weight*combTagEffCorr);
     hw1TagPartPt->Fill(ptSB,weight*combTagEffCorr);
     hw1TagPartEta->Fill(jtetaSB,weight*combTagEffCorr);
     hw1TagCent->Fill(bin/2,weight*combTagEffCorr);

     hTagXj->Fill(ptSB/pt1,weight);
     hw1TagXj->Fill(ptSB/pt1,weight*combTagEffCorr);

     totCounts3+=weight;     
     totCorr3+=weight*combTagEffCorr;
     
     
  }

 cout<<" totCorr3 / totCounts3 "<<totCorr3/totCounts3<<endl;

  hw1TagLeadPt->Scale(totCounts3/totCorr3);
  hw1TagLeadEta->Scale(totCounts3/totCorr3);
  hw1TagPartPt->Scale(totCounts3/totCorr3);
  hw1TagPartEta->Scale(totCounts3/totCorr3);
  hw1TagCent->Scale(totCounts3/totCorr3);
  hw1TagXj->Scale(totCounts3/totCorr3);
 
  hw1TagEffLeadPt->Divide(hw1TagLeadPt,hLeadPt,1.,1.,"B");  
  hw1TagEffLeadEta->Divide(hw1TagLeadEta,hLeadEta,1.,1.,"B");  
  hw1TagEffPartPt->Divide(hw1TagPartPt,hPartPt,1.,1.,"B");  
  hw1TagEffPartEta->Divide(hw1TagPartEta,hPartEta,1.,1.,"B");
  hw1TagEffPartCent->Divide(hw1TagCent,hCent,1.,1.,"B");  
  //hw1TagEffXj->Divide(hw1TagXj,hXj,1.,1.,"B");  


  if(doResid){
    
    hw1TagEffLeadPt->Fit(fLeadPt2);
    
    double totCounts4 =0;
    double totCorr4 =0;
    
    for (Long64_t i=0; i<nentries;i++) {
      nt->GetEntry(i);
      
      if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
      else if(!requireFCR && bProdCode==2) continue;
      
      if(cBin==1){
	if(bin>=20) continue;
      }
      else if(cBin==2){
	if(bin<20 || bin>=60) continue;       
      }
      else if(cBin==3){
	if(bin<60) continue;       
      }
      
      
      float pt1=jtpt1;
      if(useGenPt) pt1 = refpt1;
      
      if(pt1<100) continue;
      if(abs(refparton_flavorForB1)!=5) continue;
      
      float ptSB=jtptSB;
      if(useGenPt) ptSB = refptSB;
      
      if(ptSB<40) continue;
      if(dphiSB1<2./3.*acos(-1.)) continue;
      if(refptSB<20) continue;
	   
      if(discr_csvV1_1<0.9) continue;
      if(discr_csvV1_SB<0.9) continue;
      
      float centTagEffCorr = 1.;
      if(ppPbPb) centTagEffCorr=1./fLeadCent->Eval(bin/2.)/fPartCent->Eval(bin/2.);

      float leadTagEffPt = fLeadPt->Eval(pt1)*fLeadPt2->Eval(pt1);
      if(pt1>maxPtLead) leadTagEffPt=fLeadPt->Eval(maxPtLead)*fLeadPt2->Eval(maxPtLead);
      float leadTagEffEta = fLeadEta->Eval(jteta1);
      float leadTagEffCorr = 1./leadTagEffPt/leadTagEffEta;
      
      float partTagEffPt = fPartPt->Eval(ptSB);
      if(ptSB>maxPtPart) partTagEffPt=fPartPt->Eval(maxPtPart);
      float partTagEffEta = fPartEta->Eval(jtetaSB);
      float partTagEffCorr = 1./partTagEffPt/partTagEffEta;
      
      float combTagEffCorr = centTagEffCorr*leadTagEffCorr*partTagEffCorr;
      
      
      hw2TagLeadPt->Fill(pt1,weight*combTagEffCorr);
      hw2TagLeadEta->Fill(jteta1,weight*combTagEffCorr);
      hw2TagPartPt->Fill(ptSB,weight*combTagEffCorr);
      hw2TagPartEta->Fill(jtetaSB,weight*combTagEffCorr); 
      hw2TagCent->Fill(bin/2.,weight*combTagEffCorr); 
      
      totCounts4+=weight;     
      totCorr4+=weight*combTagEffCorr;
      
      
    }
    
    cout<<" totCorr4 / totCounts4 "<<totCorr4/totCounts4<<endl;
    cout<<" totCounts4 "<<totCounts4<<endl;
    
    hw2TagLeadPt->Scale(totCounts4/totCorr4);
    hw2TagLeadEta->Scale(totCounts4/totCorr4);
    hw2TagPartPt->Scale(totCounts4/totCorr4);
    hw2TagPartEta->Scale(totCounts4/totCorr4);
    hw2TagCent->Scale(totCounts4/totCorr4);
    
    hw2TagEffLeadPt->Divide(hw2TagLeadPt,hLeadPt,1.,1.,"B");  
    hw2TagEffLeadEta->Divide(hw2TagLeadEta,hLeadEta,1.,1.,"B");  
    hw2TagEffPartPt->Divide(hw2TagPartPt,hPartPt,1.,1.,"B");  
    hw2TagEffPartEta->Divide(hw2TagPartEta,hPartEta,1.,1.,"B");  
    hw2TagEffCent->Divide(hw2TagCent,hCent,1.,1.,"B");  



    hw2TagEffPartPt->Fit(fPartPt2);
    
    double totCounts5 =0;
    double totCorr5 =0;
    
    
    for (Long64_t i=0; i<nentries;i++) {
      nt->GetEntry(i);
      
      if(pthat<50) continue;
     if(refpt1<50) continue;
     if(requireFCR && bProdCode!=1) continue;
      else if(!requireFCR && bProdCode==2) continue;
      
      if(cBin==1){
	if(bin>=20) continue;
      }
      else if(cBin==2){
	if(bin<20 || bin>=60) continue;       
      }
      else if(cBin==3){
	if(bin<60) continue;       
      }
      
      
      float pt1=jtpt1;
      if(useGenPt) pt1 = refpt1;
      
      if(pt1<100) continue;
      if(abs(refparton_flavorForB1)!=5) continue;
      
      float ptSB=jtptSB;
      if(useGenPt) ptSB = refptSB;
      
      if(ptSB<40) continue;
      if(dphiSB1<2./3.*acos(-1.)) continue;
      
      
      if(discr_csvV1_1<0.9) continue;
      if(discr_csvV1_SB<0.9) continue;
      
      float centTagEffCorr = 1.;
      if(ppPbPb) centTagEffCorr=1./fLeadCent->Eval(bin/2.)/fPartCent->Eval(bin/2.);
      
      float leadTagEffPt = fLeadPt->Eval(pt1)*fLeadPt2->Eval(pt1);
      if(pt1>maxPtLead) leadTagEffPt=fLeadPt->Eval(maxPtLead)*fLeadPt2->Eval(maxPtLead);
      float leadTagEffEta = fLeadEta->Eval(jteta1);
      float leadTagEffCorr = 1./leadTagEffPt/leadTagEffEta;
      
      float partTagEffPt = fPartPt->Eval(ptSB)*fPartPt2->Eval(ptSB);
      if(ptSB>maxPtPart) partTagEffPt=fPartPt->Eval(maxPtPart)*fPartPt2->Eval(maxPtPart);
      float partTagEffEta = fPartEta->Eval(jtetaSB);
      float partTagEffCorr = 1./partTagEffPt/partTagEffEta;
      
      float combTagEffCorr = centTagEffCorr*leadTagEffCorr*partTagEffCorr;
      
      
      hw3TagLeadPt->Fill(pt1,weight*combTagEffCorr);
      hw3TagLeadEta->Fill(jteta1,weight*combTagEffCorr);
      hw3TagPartPt->Fill(ptSB,weight*combTagEffCorr);
      hw3TagPartEta->Fill(jtetaSB,weight*combTagEffCorr); 
      hw3TagCent->Fill(bin/2.,weight*combTagEffCorr); 
      
      totCounts5+=weight;     
      totCorr5+=weight*combTagEffCorr;
      
      
    }
    
    cout<<" totCorr5 / totCounts5 "<<totCorr5/totCounts5<<endl;
    cout<<" totCounts5 "<<totCounts5<<endl;
    
    
    hw3TagLeadPt->Scale(totCounts5/totCorr5);
    hw3TagLeadEta->Scale(totCounts5/totCorr5);
    hw3TagPartPt->Scale(totCounts5/totCorr5);
    hw3TagPartEta->Scale(totCounts5/totCorr5);
    hw3TagCent->Scale(totCounts5/totCorr5);
    
    hw3TagEffLeadPt->Divide(hw3TagLeadPt,hLeadPt,1.,1.,"B");  
    hw3TagEffLeadEta->Divide(hw3TagLeadEta,hLeadEta,1.,1.,"B");  
    hw3TagEffPartPt->Divide(hw3TagPartPt,hPartPt,1.,1.,"B");  
    hw3TagEffPartEta->Divide(hw3TagPartEta,hPartEta,1.,1.,"B");
    hw3TagEffCent->Divide(hw3TagCent,hCent,1.,1.,"B");  
  }
  
  if(doRebin){
    hTagEffLeadPt->Rebin(2);
    hTagEffLeadEta->Rebin(2);
    hTagEffPartPt->Rebin(2);
    hTagEffPartEta->Rebin(2);
    hTagEffLeadCent->Rebin(2);
    hTagEffPartCent->Rebin(2);
    
    hTagEffLeadPt->Scale(1./2.);
    hTagEffLeadEta->Scale(1./2.);
    hTagEffPartPt->Scale(1./2.);
    hTagEffPartEta->Scale(1./2.);
    hTagEffLeadCent->Scale(1./2.);
    hTagEffPartCent->Scale(1./2.);
    
    hw1TagEffLeadPt->Rebin(2);
    hw1TagEffLeadEta->Rebin(2);
    hw1TagEffPartPt->Rebin(2);
    hw1TagEffPartEta->Rebin(2);
    hw1TagEffPartCent->Rebin(2);

    hw1TagEffLeadPt->Scale(1./2.);
    hw1TagEffLeadEta->Scale(1./2.);
    hw1TagEffPartPt->Scale(1./2.);
    hw1TagEffPartEta->Scale(1./2.);
    hw1TagEffPartCent->Scale(1./2.);

    hw2TagEffLeadPt->Rebin(2);
    hw2TagEffLeadEta->Rebin(2);
    hw2TagEffPartPt->Rebin(2);
    hw2TagEffPartEta->Rebin(2);
    hw2TagEffCent->Rebin(2);

    hw2TagEffLeadPt->Scale(1./2.);
    hw2TagEffLeadEta->Scale(1./2.);
    hw2TagEffPartPt->Scale(1./2.);
    hw2TagEffPartEta->Scale(1./2.);
    hw2TagEffCent->Scale(1./2.);

    hw3TagEffLeadPt->Rebin(2);
    hw3TagEffLeadEta->Rebin(2);
    hw3TagEffPartPt->Rebin(2);
    hw3TagEffPartEta->Rebin(2);
    hw3TagEffCent->Rebin(2);
    
    hw3TagEffLeadPt->Scale(1./2.);
    hw3TagEffLeadEta->Scale(1./2.);
    hw3TagEffPartPt->Scale(1./2.);
    hw3TagEffPartEta->Scale(1./2.);
    hw3TagEffCent->Scale(1./2.);

  }

  hw1TagEffLeadPt->SetMarkerColor(2);
  hw1TagEffLeadEta->SetMarkerColor(2);
  hw1TagEffPartPt->SetMarkerColor(2);
  hw1TagEffPartEta->SetMarkerColor(2);
  hw1TagEffPartCent->SetMarkerColor(2);

  hw1TagEffLeadPt->SetLineColor(2);
  hw1TagEffLeadEta->SetLineColor(2);
  hw1TagEffPartPt->SetLineColor(2);
  hw1TagEffPartEta->SetLineColor(2);
  hw1TagEffPartCent->SetLineColor(2);

  
  hw2TagEffLeadPt->SetMarkerColor(3);
  hw2TagEffLeadEta->SetMarkerColor(3);
  hw2TagEffPartPt->SetMarkerColor(3);
  hw2TagEffPartEta->SetMarkerColor(3);
  hw2TagEffCent->SetMarkerColor(3);

  hw2TagEffLeadPt->SetLineColor(3);
  hw2TagEffLeadEta->SetLineColor(3);
  hw2TagEffPartPt->SetLineColor(3);
  hw2TagEffPartEta->SetLineColor(3);
  hw2TagEffCent->SetLineColor(3);
  
  hw3TagEffLeadPt->SetMarkerColor(4);
  hw3TagEffLeadEta->SetMarkerColor(4);
  hw3TagEffPartPt->SetMarkerColor(4);
  hw3TagEffPartEta->SetMarkerColor(4);
  hw3TagEffCent->SetMarkerColor(4);

  hw3TagEffLeadPt->SetLineColor(4);
  hw3TagEffLeadEta->SetLineColor(4);
  hw3TagEffPartPt->SetLineColor(4);
  hw3TagEffPartEta->SetLineColor(4);
  hw3TagEffCent->SetLineColor(4);

  TLegend *leg=new TLegend(0.6,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hTagEffLeadPt,"Uncorrected","p");
  leg->AddEntry(hw1TagEffLeadPt,"Corrected","p");

  TLegend *leg2=new TLegend(0.2,0.7,0.5,0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hTagEffLeadPt,"Leading","p");
  leg2->AddEntry(hTagEffPartPt,"Subleading","p");
  
  
  TCanvas *c0=new TCanvas("c0","c0",600,600);
  hTagEffLeadCent->SetMinimum(0.);
  hTagEffLeadCent->SetMaximum(1.);
  hTagEffLeadCent->Draw();
  //  hTagEffPartCent->SetMarkerStyle(4);
  //  hTagEffPartCent->Draw("same");
  hw1TagEffPartCent->Draw("same");
  if(doResid){
    hw2TagEffCent->Draw("same");
    hw3TagEffCent->Draw("same");
  }
  leg->Draw();
  TF1 *line0 = new TF1("line0","pol0",0,100);
  line0->SetLineColor(2);
  line0->SetLineWidth(1);
  line0->SetLineStyle(2);
  hw1TagEffPartCent->Fit(line0);


  TCanvas *c00=new TCanvas("c00","c00",600,600);
  hTagEffPartCent->SetMinimum(0.);
  hTagEffPartCent->SetMaximum(1.);
  hTagEffPartCent->Draw();
  hw1TagEffPartCent->Draw("same");
  if(doResid){
    hw2TagEffCent->Draw("same");
    hw3TagEffCent->Draw("same");
  }
  leg->Draw();
  TF1 *line00 = new TF1("line0","pol0",0,100);
  line00->SetLineColor(2);
  line00->SetLineWidth(1);
  line00->SetLineStyle(2);
  hw1TagEffPartCent->Fit(line0);



  
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  hTagEffLeadPt->SetMinimum(0.);
  hTagEffLeadPt->SetMaximum(1.);
  hTagEffLeadPt->Draw();
  hw1TagEffLeadPt->Draw("same");
  if(doResid){
    hw2TagEffLeadPt->Draw("same");
    hw3TagEffLeadPt->Draw("same");
  }
  leg->Draw();
  TF1 *line1 = new TF1("line1","pol0",100,maxPtLead);
  line1->SetLineColor(2);
  line1->SetLineWidth(1);
  line1->SetLineStyle(2);
  hw1TagEffLeadPt->Fit(line1);

  
  TCanvas *c2=new TCanvas("c2","c2",600,600);
  hTagEffLeadEta->SetMinimum(0.);
  hTagEffLeadEta->SetMaximum(1.);
  hTagEffLeadEta->Draw();
  hw1TagEffLeadEta->Draw("same");
    if(doResid){
      hw2TagEffLeadEta->Draw("same");
      hw3TagEffLeadEta->Draw("same");
    }
    leg->Draw();
  TF1 *line2 = new TF1("line2","pol0",-2,2);
  line2->SetLineColor(2);
  line2->SetLineWidth(1);
  line2->SetLineStyle(2);
  hw1TagEffLeadEta->Fit(line2);

    TCanvas *c3=new TCanvas("c3","c3",600,600);
  hTagEffPartPt->SetMinimum(0.);
  hTagEffPartPt->SetMaximum(1.);
  hTagEffPartPt->Draw();
  hw1TagEffPartPt->Draw("same");
  if(doResid){
    hw2TagEffPartPt->Draw("same");
    hw3TagEffPartPt->Draw("same");
  }
  leg->Draw();
  TF1 *line3 = new TF1("line3","pol0",40,maxPtPart);
  line3->SetLineColor(2);
  line3->SetLineWidth(1);
  line3->SetLineStyle(2);
  hw1TagEffPartPt->Fit(line3);
    
    TCanvas *c4=new TCanvas("c4","c4",600,600);
    hTagEffPartEta->SetMinimum(0.);
    hTagEffPartEta->SetMaximum(1.);
    hTagEffPartEta->Draw();
    hw1TagEffPartEta->Draw("same");
    if(doResid){
      hw2TagEffPartEta->Draw("same");
      hw3TagEffPartEta->Draw("same");
    }
    leg->Draw();
    TF1 *line4 = new TF1("line4","pol0",-2,2);
  line4->SetLineColor(2);
  line4->SetLineWidth(1);
  line4->SetLineStyle(2);
  hw1TagEffPartEta->Fit(line4);


  TCanvas *c5=new TCanvas("c5","c5",600,600);

  hXj->SetMarkerStyle(24);
  hTagXj->SetMarkerStyle(25);
  hw1TagXj->SetMarkerStyle(26);
  hXj->SetMarkerColor(1);
  hTagXj->SetMarkerColor(2);
  hw1TagXj->SetMarkerColor(3);
  
  hXj->Scale(1./hXj->Integral());
  hTagXj->Scale(1./hTagXj->Integral());
  hw1TagXj->Scale(1./hw1TagXj->Integral());
  hXj->Draw();
   hTagXj->Draw("same");
   hw1TagXj->Draw("same");
   leg->Draw();
  
   cout<<" Mean xJ = "<<hXj->GetMean()<<endl;
   cout<<" Mean xJ tagged = "<<hTagXj->GetMean()<<endl;
   cout<<" Mean xJ corrected = "<<hw1TagXj->GetMean()<<endl;

  
    TFile *fout = NULL;
  if(ppPbPb){    
    if(requireFCR) fout=new TFile(Form("offTagEff_PbPb_fcr_cBin_%i_v3.root",cBin),"recreate");
    else fout=new TFile(Form("offTagEff_PbPb_all_cBin%i_v3.root",cBin),"recreate");
  }
  else{
    if(requireFCR) fout=new TFile("offTagEff_pp_fcr_v3.root","recreate");
    else fout=new TFile("offTagEff_pp_all_v3.root","recreate");
  }

  fout->cd();
  /*
  hLeadPt->Write();
  hLeadEta->Write();
  hTagLeadPt->Write();
  hTagLeadEta->Write();
  hTagEffLeadPt->Write();
  hTagEffLeadEta->Write();

  hPartPt->Write();
  hPartEta->Write();
  hTagPartPt->Write();
  hTagPartEta->Write();
  hTagEffPartPt->Write();
  hTagEffPartEta->Write();
  */

  hCent->Write();
  hw1TagEffPartCent->Write();
  hw1TagCent->Write();
  hLeadPt->Write();
  hw1TagLeadPt->Write();
  hw1TagEffLeadPt->Write();


  if(ppPbPb){
    fLeadCent->Write();
    fPartCent->Write();
  }
  fLeadPt->Write();
  fLeadEta->Write();
  fPartPt->Write();
  fPartEta->Write();
  if(doResid){
    fLeadPt2->Write();
    fPartPt2->Write();  
  }
  fout->Close();

  if(savePlots){
    
    
    if(ppPbPb){
      c0->SaveAs(Form("savedPlots/tagEffVsLeadCent_cBin_%i.pdf",cBin));
      c0->SaveAs(Form("savedPlots/tagEffVsLeadCent_cBin_%i.gif",cBin));
      c0->SaveAs(Form("savedPlots/tagEffVsLeadCent_cBin_%i.C",cBin));
      c00->SaveAs(Form("savedPlots/tagEffVsPartCent_cBin_%i.pdf",cBin));
      c00->SaveAs(Form("savedPlots/tagEffVsPartCent_cBin_%i.gif",cBin));
      c00->SaveAs(Form("savedPlots/tagEffVsPartCent_cBin_%i.C",cBin));
      c1->SaveAs(Form("savedPlots/tagEffVsLeadPt_cBin_%i.pdf",cBin));
      c1->SaveAs(Form("savedPlots/tagEffVsLeadPt_cBin_%i.gif",cBin));
      c1->SaveAs(Form("savedPlots/tagEffVsLeadPt_cBin_%i.C",cBin));
      c2->SaveAs(Form("savedPlots/tagEffVsLeadEta_cBin_%i.pdf",cBin));
      c2->SaveAs(Form("savedPlots/tagEffVsLeadEta_cBin_%i.gif",cBin));
      c2->SaveAs(Form("savedPlots/tagEffVsLeadEta_cBin_%i.C",cBin));
      c3->SaveAs(Form("savedPlots/tagEffVsPartPt_cBin_%i.pdf",cBin));
      c3->SaveAs(Form("savedPlots/tagEffVsPartPt_cBin_%i.gif",cBin));
      c3->SaveAs(Form("savedPlots/tagEffVsPartPt_cBin_%i.C",cBin));
      c4->SaveAs(Form("savedPlots/tagEffVsPartEta_cBin_%i.pdf",cBin));
      c4->SaveAs(Form("savedPlots/tagEffVsPartEta_cBin_%i.gif",cBin));
      c4->SaveAs(Form("savedPlots/tagEffVsPartEta_cBin_%i.C",cBin));
    }
    else{
      c1->SaveAs("savedPlots/tagEffVsLeadPt.pdf");
      c1->SaveAs("savedPlots/tagEffVsLeadPt.gif");
      c1->SaveAs("savedPlots/tagEffVsLeadPt.C");
      c2->SaveAs("savedPlots/tagEffVsLeadEta.pdf");
      c2->SaveAs("savedPlots/tagEffVsLeadEta.gif");
      c2->SaveAs("savedPlots/tagEffVsLeadEta.C");
      c3->SaveAs("savedPlots/tagEffVsPartPt.pdf");
      c3->SaveAs("savedPlots/tagEffVsPartPt.gif");
      c3->SaveAs("savedPlots/tagEffVsPartPt.C");
      c4->SaveAs("savedPlots/tagEffVsPartEta.pdf");
      c4->SaveAs("savedPlots/tagEffVsPartEta.gif");
      c4->SaveAs("savedPlots/tagEffVsPartEta.C");
      
    }
    
  }


}
