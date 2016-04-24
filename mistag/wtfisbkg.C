#include "../helpers/plotting.h"
#include "../helpers/config.h"

void RenameBinLabelsX(TH1 *h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    h->GetXaxis()->SetBinLabel(i,Form("%d",i-1));
}

void RenameBinLabelsX(TH1 *h, vector<TString> labels)
{
  if (h->GetNbinsX()!=labels.size()) cout<<" wrong number of labels/bins"<<endl;

  for (int i=1;i<=h->GetNbinsX();i++)
    h->GetXaxis()->SetBinLabel(i,labels[i-1]);
}

void RenameBinLabelsY(TH1 *h)
{
  for (int i=1;i<=h->GetNbinsY();i++)
    h->GetYaxis()->SetBinLabel(i,Form("%d",i-1));
}

void RenameBinLabelsY(TH1 *h, vector<TString> labels)
{
    if (h->GetNbinsY()!=labels.size()) cout<<" wrong number of labels/bins"<<endl;
  for (int i=1;i<=h->GetNbinsY();i++)
    h->GetYaxis()->SetBinLabel(i,labels[i-1]);
}

void wtfisbkg()
{
  auto fmcpp = new TFile(config.getFileName_djt("mcppbfa"));
  auto nt = (TTree *)fmcpp->Get("nt");

  buildh(4,0,4);

  auto bprodcode45 = geth("bprodcode45","C-B;bProdCode");
  auto bprodcode54 = geth("bprodcode54","B-C;bProdCode");
  nt->Project("bprodcode45","bProdCode","weight*(abs(refparton_flavorForB1)==4 && abs(refparton_flavorForBSL)==5 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("bprodcode54","bProdCode","weight*(abs(refparton_flavorForB1)==5 && abs(refparton_flavorForBSL)==4 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  
  buildh(7,0,7,7,0,7);
  auto flavproc54_0 = geth2d("flavproc54_0","B-C, GSP;flavorProcess1 (B);flavorProcessSL (C)");
  auto flavproc54_1 = geth2d("flavproc54_1","B-C, FCR;flavorProcess1 (B);flavorProcessSL (C)");
  auto flavproc54_2 = geth2d("flavproc54_2","B-C, FEX;flavorProcess1 (B);flavorProcessSL (C)");
  auto flavproc45_0 = geth2d("flavproc45_0","C-B, GSP;flavorProcess1 (C);flavorProcessSL (B)");
  auto flavproc45_1 = geth2d("flavproc45_1","C-B, FCR;flavorProcess1 (C);flavorProcessSL (B)");
  auto flavproc45_2 = geth2d("flavproc45_2","C-B, FEX;flavorProcess1 (C);flavorProcessSL (B)");

  nt->Project("flavproc54_0", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==0 && abs(refparton_flavorForB1)==5 && abs(refparton_flavorForBSL)==4 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("flavproc54_1", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==1 && abs(refparton_flavorForB1)==5 && abs(refparton_flavorForBSL)==4 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("flavproc54_2", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==2 && abs(refparton_flavorForB1)==5 && abs(refparton_flavorForBSL)==4 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("flavproc45_0", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==0 && abs(refparton_flavorForB1)==4 && abs(refparton_flavorForBSL)==5 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("flavproc45_1", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==1 && abs(refparton_flavorForB1)==4 && abs(refparton_flavorForBSL)==5 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");
  nt->Project("flavproc45_2", "refparton_flavorProcessSL:refparton_flavorProcess1","weight*(bProdCode==2 && abs(refparton_flavorForB1)==4 && abs(refparton_flavorForBSL)==5 && jtpt1>100 && refpt1>50 && pthat>50 && jtptSL>40 && dphiSL1>2.1 && discr_csvV1_1>0.9)");


  RenameBinLabelsX(bprodcode45,{"GSP","FCR","FEX","FEX2"});
  RenameBinLabelsX(bprodcode54,{"GSP","FCR","FEX","FEX2"});
  
  RenameBinLabelsX(flavproc54_0,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsX(flavproc54_1,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsX(flavproc54_2,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsX(flavproc45_0,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsX(flavproc45_1,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsX(flavproc45_2,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc54_0,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc54_1,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc54_2,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc45_0,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc45_1,{"","prim","","interm","interm","final","final"});
  RenameBinLabelsY(flavproc45_2,{"","prim","","interm","interm","final","final"});


  aktstring = "";

  Draw({bprodcode45,bprodcode54});

  Draw({flavproc54_0},"colz");
  Draw({flavproc54_1},"colz");
  Draw({flavproc54_2},"colz");
  Draw({flavproc45_0},"colz");
  Draw({flavproc45_1},"colz");
  Draw({flavproc45_2},"colz");


}
