bool applyCorrection = true;

#include "../helpers/config.h"

class tageffcorr 
{
public:
  vector<TF1 *> fPbPbLeadCent,fPbPbPartCent,fPbPbLeadPt,fPbPbLeadEta,fPbPbPartPt,fPbPbPartEta,fPbPbLeadPt2,fPbPbPartPt2;
  TF1 *fLeadPt,*fLeadEta,*fPartPt,*fPartEta,*fLeadPt2,*fPartPt2;

  TString tagcorpath="";



  tageffcorr(float csv)
  {
    // map<float,TString> csvfolder = {{0.9,"tagging_jtsignal2v4_0p90"},
    //                                 {0.85,"tagging_jtsignal2v4_0p85"},
    //                                 {0.95,"tagging_jtsignal2v4_0p95"}};
 map<float,TString> csvfolder = {
    {0.00,"tagging_jtsignal2v4_0p00"},
    {0.10,"tagging_jtsignal2v4_0p10"},
    {0.20,"tagging_jtsignal2v4_0p20"},
    {0.30,"tagging_jtsignal2v4_0p30"},
    {0.40,"tagging_jtsignal2v4_0p40"},
    {0.50,"tagging_jtsignal2v4_0p50"},
    {0.60,"tagging_jtsignal2v4_0p60"},
    {0.70,"tagging_jtsignal2v4_0p70"},
    {0.80,"tagging_jtsignal2v4_0p80"},
    {0.85,"tagging_jtsignal2v4_0p85"},
    {0.90,"tagging_jtsignal2v4_0p90"},
    {0.95,"tagging_jtsignal2v4_0p95"},
  };
                                    

    if (csvfolder.find(csv)==csvfolder.end()) {
        cout<<"Unknown correction for csv="<<csv<<endl;
        return;
    }

    tagcorpath = config.tagcorfolder+"/"+csvfolder[csv];

    cout<<"Using csv = "<<csv<<", tagging corrections at "<<tagcorpath<<endl;

    loadTagEffCorrections();

  }

  void loadTagEffCorrections();
  float pp(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta);
  float PbPb(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin);

  float pp1(float leadingjetpt, float leadingjeteta);
  float pp2(float partnerjetpt, float partnerjeteta);
  float PbPb1(float leadingjetpt, float leadingjeteta, float bin);
  float PbPb2(float partnerjetpt, float partnerjeteta, float bin);


};


void tageffcorr::loadTagEffCorrections()
{
  fPbPbLeadCent.resize(3);
  fPbPbPartCent.resize(3);
  fPbPbLeadPt.resize(3);
  fPbPbLeadEta.resize(3);
  fPbPbPartPt.resize(3);
  fPbPbPartEta.resize(3);


  TFile *fppFits       = new TFile(tagcorpath+"/offTagEff_pp_all_v4.root");
  TFile *fPbPbFitsbin1 = new TFile(tagcorpath+"/offTagEff_PbPb_all_cBin1_v4.root");
  TFile *fPbPbFitsbin2 = new TFile(tagcorpath+"/offTagEff_PbPb_all_cBin2_v4.root");
  TFile *fPbPbFitsbin3 = new TFile(tagcorpath+"/offTagEff_PbPb_all_cBin3_v4.root");

  fLeadPt =  (TF1 *)fppFits->Get("piecewisePt1");
  fLeadEta = (TF1 *)fppFits->Get("fLeadEta");
  fPartPt =  (TF1 *)fppFits->Get("fPartPt");
  fPartEta = (TF1 *)fppFits->Get("fPartEta");

  fPbPbLeadCent[0]  = (TF1 *)fPbPbFitsbin1->Get("fLeadCent");
  fPbPbPartCent[0]  = (TF1 *)fPbPbFitsbin1->Get("fPartCent");
  fPbPbLeadPt[0]  = (TF1 *)fPbPbFitsbin1->Get("piecewisePt1");
  fPbPbLeadEta[0] = (TF1 *)fPbPbFitsbin1->Get("fLeadEta");
  fPbPbPartPt[0]  = (TF1 *)fPbPbFitsbin1->Get("fPartPt");
  fPbPbPartEta[0] = (TF1 *)fPbPbFitsbin1->Get("fPartEta");
  
  fPbPbLeadCent[1]  = (TF1 *)fPbPbFitsbin2->Get("fLeadCent");
  fPbPbPartCent[1]  = (TF1 *)fPbPbFitsbin2->Get("fPartCent");
  fPbPbLeadPt[1]  = (TF1 *)fPbPbFitsbin2->Get("piecewisePt1");
  fPbPbLeadEta[1] = (TF1 *)fPbPbFitsbin2->Get("fLeadEta");
  fPbPbPartPt[1]  = (TF1 *)fPbPbFitsbin2->Get("fPartPt");
  fPbPbPartEta[1] = (TF1 *)fPbPbFitsbin2->Get("fPartEta");
  
  fPbPbLeadCent[2]  = (TF1 *)fPbPbFitsbin3->Get("fLeadCent");
  fPbPbPartCent[2]  = (TF1 *)fPbPbFitsbin3->Get("fPartCent");
  fPbPbLeadPt[2]  = (TF1 *)fPbPbFitsbin3->Get("piecewisePt1");
  fPbPbLeadEta[2] = (TF1 *)fPbPbFitsbin3->Get("fLeadEta");
  fPbPbPartPt[2]  = (TF1 *)fPbPbFitsbin3->Get("fPartPt");
  fPbPbPartEta[2] = (TF1 *)fPbPbFitsbin3->Get("fPartEta");


}



float tageffcorr::pp(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta)
{
  return pp1(leadingjetpt,leadingjeteta)*pp2(partnerjetpt,partnerjeteta);
}

float tageffcorr::pp1(float leadingjetpt, float leadingjeteta)
{

  if (!applyCorrection) return 1;
  if (leadingjetpt>200) leadingjetpt = 200;
  if (leadingjetpt<100) leadingjetpt = 100;
  
  return 1/(fLeadPt->Eval(leadingjetpt)*fLeadEta->Eval(leadingjeteta));
}

float tageffcorr::pp2(float partnerjetpt, float partnerjeteta)
{

  if (!applyCorrection) return 1;
  if (partnerjetpt>140) partnerjetpt = 140;
  if (partnerjetpt<40) partnerjetpt = 40;

  return 1/(fPartPt->Eval(partnerjetpt)*fPartEta->Eval(partnerjeteta));
}


float tageffcorr::PbPb(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin)
{

//   if (!applyCorrection) return 1;

//   int b = getbinindex(bin);

//   if (leadingjetpt>200) leadingjetpt = 200;
//   if (partnerjetpt>140) partnerjetpt = 140;

//   double corr = 1/fPbPbLeadCent[b]->Eval(bin/2)/fPbPbPartCent[b]->Eval(bin/2)//1/fPbPbCent[b]->Eval(bin/2)
//                                          /fPbPbLeadPt[b]->Eval(leadingjetpt)/ fPbPbLeadEta[b]->Eval(leadingjeteta)
//                                          /fPbPbPartPt[b]->Eval(partnerjetpt)
//                                          /fPbPbPartEta[b]->Eval(partnerjeteta);

// //                 /fPbPbLeadPt2[b]->Eval(leadingjetpt)/fPbPbPartPt2[b]->Eval(partnerjetpt);
  return PbPb1(leadingjetpt,leadingjeteta,bin)*PbPb2(partnerjetpt,partnerjeteta,bin);
}


float tageffcorr::PbPb1(float leadingjetpt, float leadingjeteta, float bin)
{
  if (!applyCorrection) return 1;
  int b = getbinindex(bin);
  if (leadingjetpt>200) leadingjetpt = 200;

  return 1/(fPbPbLeadCent[b]->Eval(bin/2)*fPbPbLeadPt[b]->Eval(leadingjetpt)*fPbPbLeadEta[b]->Eval(leadingjeteta));
}


float tageffcorr::PbPb2(float partnerjetpt, float partnerjeteta, float bin)
{
  if (!applyCorrection) return 1;
  int b = getbinindex(bin);
  if (partnerjetpt>140) partnerjetpt = 140;

  return 1/(fPbPbPartCent[b]->Eval(bin/2)*fPbPbPartPt[b]->Eval(partnerjetpt)*fPbPbPartEta[b]->Eval(partnerjeteta));

}



tageffcorr *tagcorj1=NULL, *tagcorj2=NULL;

void inittageffcorr(float csv1, float csv2)
{
  csvcut1 = csv1;
  csvcut2 = csv2;

  tagcorj1 = new tageffcorr(csv1);
  tagcorj2 = new tageffcorr(csv2);
}

float tageffcorrectionPbPb(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin)
{
  return tagcorj1->PbPb1(leadingjetpt,leadingjeteta,bin)*tagcorj2->PbPb2(partnerjetpt,partnerjeteta,bin);
}

float tageffcorrectionpp(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta)
{
  return tagcorj1->pp1(leadingjetpt,leadingjeteta)*tagcorj2->pp2(partnerjetpt,partnerjeteta);
}
