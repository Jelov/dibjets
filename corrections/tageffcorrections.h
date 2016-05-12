bool applyCorrection = true;

TF1 *fppEta, *fppPt;
//TF1 *fPbPbCent, *fPbPbEta, *fPbPbPt;
TF1 *fPbPbCent0, *fPbPbEta0, *fPbPbPt0;
TF1 *fPbPbCent1, *fPbPbEta1, *fPbPbPt1;
TF1 *fPbPbCent2, *fPbPbEta2, *fPbPbPt2;

vector<TF1 *> //fPbPbCent(3),
              fPbPbLeadCent(3),fPbPbPartCent(3),fPbPbLeadPt(3),fPbPbLeadEta(3),fPbPbPartPt(3),fPbPbPartEta(3),fPbPbLeadPt2(3),fPbPbPartPt2(3);

TF1 *fLeadPt,*fLeadEta,*fPartPt,*fPartEta,*fLeadPt2,*fPartPt2;

void loadTagEffCorrections()
{
  TFile *fppFits = new TFile("../../ntuples/offTagEff_pp_all.root");//offlineTagEff_PP.root");
   
   // fppEta = (TF1*) fppFits->Get("fitEta");
   // fppPt = (TF1*) fppFits->Get("fitPt");

  fLeadPt =  (TF1 *)fppFits->Get("piecewisePt1");
  fLeadEta = (TF1 *)fppFits->Get("fLeadEta");
  fPartPt =  (TF1 *)fppFits->Get("fPartPt");
  fPartEta = (TF1 *)fppFits->Get("fPartEta");
  // fLeadPt2 = (TF1 *)fppFits->Get("fLeadPt2");
  // fPartPt2 = (TF1 *)fppFits->Get("fPartPt2");


  // TFile *fPbPbFits = new TFile("../../ntuples/offlineTagEff_PbPb_take5.root");//offlineTagEff_PbPb.root");
  // fPbPbCent = (TF1*) fPbPbFits->Get("fitCent");
  // fPbPbEta = (TF1*) fPbPbFits->Get("fitEta");
  // fPbPbPt = (TF1*) fPbPbFits->Get("fitPt");

  TFile *fPbPbFitsbin1 = new TFile("../../ntuples/offTagEff_PbPb_all_cBin1_v2.root");
  TFile *fPbPbFitsbin2 = new TFile("../../ntuples/offTagEff_PbPb_all_cBin2_v2.root");
  TFile *fPbPbFitsbin3 = new TFile("../../ntuples/offTagEff_PbPb_all_cBin3_v2.root");

// fPbPbCent[0]  = (TF1 *)fPbPbFitsbin1->Get("fCent");
fPbPbLeadCent[0]  = (TF1 *)fPbPbFitsbin1->Get("fLeadCent");
fPbPbPartCent[0]  = (TF1 *)fPbPbFitsbin1->Get("fPartCent");
fPbPbLeadPt[0]  = (TF1 *)fPbPbFitsbin1->Get("piecewisePt1");
fPbPbLeadEta[0] = (TF1 *)fPbPbFitsbin1->Get("fLeadEta");
fPbPbPartPt[0]  = (TF1 *)fPbPbFitsbin1->Get("fPartPt");
fPbPbPartEta[0] = (TF1 *)fPbPbFitsbin1->Get("fPartEta");
// fPbPbLeadPt2[0] = (TF1 *)fPbPbFitsbin1->Get("fLeadPt2");
// fPbPbPartPt2[0] = (TF1 *)fPbPbFitsbin1->Get("fPartPt2");

// fPbPbCent[1]  = (TF1 *)fPbPbFitsbin2->Get("fCent");
fPbPbLeadCent[1]  = (TF1 *)fPbPbFitsbin2->Get("fLeadCent");
fPbPbPartCent[1]  = (TF1 *)fPbPbFitsbin2->Get("fPartCent");
fPbPbLeadPt[1]  = (TF1 *)fPbPbFitsbin2->Get("piecewisePt1");
fPbPbLeadEta[1] = (TF1 *)fPbPbFitsbin2->Get("fLeadEta");
fPbPbPartPt[1]  = (TF1 *)fPbPbFitsbin2->Get("fPartPt");
fPbPbPartEta[1] = (TF1 *)fPbPbFitsbin2->Get("fPartEta");
// fPbPbLeadPt2[1] = (TF1 *)fPbPbFitsbin2->Get("fLeadPt2");
// fPbPbPartPt2[1] = (TF1 *)fPbPbFitsbin2->Get("fPartPt2");

// fPbPbCent[2]  = (TF1 *)fPbPbFitsbin3->Get("fCent");
fPbPbLeadCent[2]  = (TF1 *)fPbPbFitsbin3->Get("fLeadCent");
fPbPbPartCent[2]  = (TF1 *)fPbPbFitsbin3->Get("fPartCent");
fPbPbLeadPt[2]  = (TF1 *)fPbPbFitsbin3->Get("piecewisePt1");
fPbPbLeadEta[2] = (TF1 *)fPbPbFitsbin3->Get("fLeadEta");
fPbPbPartPt[2]  = (TF1 *)fPbPbFitsbin3->Get("fPartPt");
fPbPbPartEta[2] = (TF1 *)fPbPbFitsbin3->Get("fPartEta");
// fPbPbLeadPt2[2] = (TF1 *)fPbPbFitsbin3->Get("fLeadPt2");
// fPbPbPartPt2[2] = (TF1 *)fPbPbFitsbin3->Get("fPartPt2");

//
/*   TFile *fPbPbFits;
   
  fPbPbFits = new TFile("/data_CMS/cms/mnguyen/bJet2015/tagEffCorr/offlineTagEff_PbPb_bin0.root");
  fPbPbCent0 = (TF1*) fPbPbFits->Get("fitCent");
  fPbPbEta0 = (TF1*) fPbPbFits->Get("fitEta");
  fPbPbPt0 = (TF1*) fPbPbFits->Get("fitPt");
  fPbPbFits->Close();

  fPbPbFits = new TFile("/data_CMS/cms/mnguyen/bJet2015/tagEffCorr/offlineTagEff_PbPb_bin1.root");
  fPbPbCent1 = (TF1*) fPbPbFits->Get("fitCent");
  fPbPbEta1 = (TF1*) fPbPbFits->Get("fitEta");
  fPbPbPt1 = (TF1*) fPbPbFits->Get("fitPt");
  fPbPbFits->Close();

  fPbPbFits = new TFile("/data_CMS/cms/mnguyen/bJet2015/tagEffCorr/offlineTagEff_PbPb_bin2.root");
  fPbPbCent2 = (TF1*) fPbPbFits->Get("fitCent");
  fPbPbEta2 = (TF1*) fPbPbFits->Get("fitEta");
  fPbPbPt2 = (TF1*) fPbPbFits->Get("fitPt");
  fPbPbFits->Close();
   */
}

// float getppcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta)
// {

//   if (!applyCorrection) return 1;

//   double leadingPtTagEffCorr = (leadingjetpt > 200) ? fppPt->Eval(200):fppPt->Eval(leadingjetpt);
//   double leadingTagEffCorr = 1./fppEta->Eval(leadingjeteta)/leadingPtTagEffCorr;


//   double parterPtTagEffCorr = (partnerjetpt > 200) ? fppPt->Eval(200):fppPt->Eval(partnerjetpt);
//   double partnerTagEffCorr = 1./fppEta->Eval(partnerjeteta)/parterPtTagEffCorr;
     
//   // numerical factor will just give an average correction of 1, can omit this
//   double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr/8.98489;
//   return combinedTagEffCorr;
// }

float getppcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta)
{

  if (!applyCorrection) return 1;

  if (leadingjetpt>200) leadingjetpt = 200;
  if (partnerjetpt>140) partnerjetpt = 140;

  double corr = 1/fLeadPt->Eval(leadingjetpt)/fLeadEta->Eval(leadingjeteta)
                 /fPartPt->Eval(partnerjetpt)/fPartEta->Eval(partnerjeteta);
//                 /fLeadPt2->Eval(leadingjetpt)/fPartPt2->Eval(partnerjetpt);

  return corr;
}

float getPbPbcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin)
{

  if (!applyCorrection) return 1;

  int b = getbinindex(bin);

  if (leadingjetpt>200) leadingjetpt = 200;
  if (partnerjetpt>140) partnerjetpt = 140;

  double corr = 1/fPbPbLeadCent[b]->Eval(bin/2)/fPbPbPartCent[b]->Eval(bin/2)//1/fPbPbCent[b]->Eval(bin/2)
                                         /fPbPbLeadPt[b]->Eval(leadingjetpt)/ fPbPbLeadEta[b]->Eval(leadingjeteta)
                                         /fPbPbPartPt[b]->Eval(partnerjetpt)
                                         /fPbPbPartEta[b]->Eval(partnerjeteta);

//                 /fPbPbLeadPt2[b]->Eval(leadingjetpt)/fPbPbPartPt2[b]->Eval(partnerjetpt);
  return corr;
}

// float getPbPbcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin)
// {

//   if (!applyCorrection) return 1;
//   /*
//   if (bin<=20) {
//     fPbPbCent = fPbPbCent0;
//     fPbPbEta = fPbPbEta0;
//     fPbPbPt = fPbPbPt0;
//   } else if (bin<=60) {
//     fPbPbCent = fPbPbCent1;
//     fPbPbEta = fPbPbEta1;
//     fPbPbPt = fPbPbPt1;
//   } else {
//     fPbPbCent = fPbPbCent2;
//     fPbPbEta = fPbPbEta2;
//     fPbPbPt = fPbPbPt2;
//   }
//   */
//   double centTagEffCorr = 1./fPbPbCent->Eval(bin/2.)/fPbPbCent->Eval(bin/2.);

//   double leadingPtTagEff = (leadingjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(leadingjetpt);
//   double leadingTagEffCorr = 1./fPbPbEta->Eval(leadingjeteta)/leadingPtTagEff;


//   double parterPtTagEff = (partnerjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(partnerjetpt);
//   double partnerTagEffCorr = 1./fPbPbEta->Eval(partnerjeteta)/parterPtTagEff;

//   double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr*centTagEffCorr/488.652;
//   return combinedTagEffCorr;
// }





// float getPbPbcorrectionPt(float leadingjetpt, float partnerjetpt)
// {return 1;

//   if (!applyCorrection) return 1;

//   double leadingPtTagEff = (leadingjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(leadingjetpt);
//   double leadingTagEffCorr = 1./leadingPtTagEff;

//   double parterPtTagEff = (partnerjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(partnerjetpt);
//   double partnerTagEffCorr = 1./parterPtTagEff;

//   double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr;

//   return combinedTagEffCorr;
// }


// float getPbPbcorrectionEta(float leadingjeteta, float partnerjeteta)
// {return 1;

//   if (!applyCorrection) return 1;

//   double leadingTagEffCorr = 1./fPbPbEta->Eval(leadingjeteta);

//   double partnerTagEffCorr = 1./fPbPbEta->Eval(partnerjeteta);

//   double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr;
//   return combinedTagEffCorr;
// }



// float getPbPbcorrectionBin(float bin)
// {return 1;

//   if (!applyCorrection) return 1;

//   double centTagEffCorr = 1./fPbPbCent->Eval(bin/2.)/fPbPbCent->Eval(bin/2.);

//   double combinedTagEffCorr =centTagEffCorr;
//   return combinedTagEffCorr;
// }

