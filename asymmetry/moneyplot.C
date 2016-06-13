#include "../helpers/looptuple.h"
#include "../helpers/plotting.h"
#include "../helpers/physics.h"
#include "./CMS_lumi.C"
#include "TExec.h"

TH1F *hmcincsys, *hdtincsys, *hmcbjtsys, *hdtbjtsys;

float addIn2(float v1, float v2)
{
  return sqrt(v1*v1+v2*v2);
}

void fitPaleLine(vector<TH1F *> v)
{
  for (auto h:v) {
    h->Fit("pol1");
    auto f= h->GetFunction("pol1");
    f->SetLineColor(h->GetLineColor());
    f->SetLineWidth(1);
    f->SetLineStyle(2);
  }
}

void moneyplot(TString name="")
{

  // Lumi stuff

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  //int iPos=11; // : top-left, left-aligned
  int  iPos=33;// : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)


  macro m("moneyplot_"+name);

  auto res = ReadFromFile("results_"+name+"/results.root");

  int nbins = binnames.size()+1;
  seth(nbins,0,nbins);
  auto hmcinc       = geth("hmcinc","MC Inclusive;;#LTx_{J}#GT");
  auto hmcsig       = geth("hmcsig","MC Inclusive Signal;;#LTx_{J}#GT");
  auto hdtinc       = geth("hdtinc","Data Inclusive;;#LTx_{J}#GT");
  auto hmcbjt       = geth("hmcbjt","MC b-jets;;#LTx_{J}#GT");
  auto hmcbSB       = geth("hmcbSB","MC b-jets Signal;;#LTx_{J}#GT");
  auto hdtbjt       = geth("hdtbjt","Data b-jets;;#LTx_{J}#GT");
  auto hdtb12       = geth("hdtb12","Data b-jets 12;;#LTx_{J}#GT");
  auto hmcb12       = geth("hmcb12","MC b-jets 12;;#LTx_{J}#GT");
  auto hmcb12Signal = geth("hmcb12Signal","MC b-jets 12 Signal;;#LTx_{J}#GT");

  
  vector<TString> labels = {"pp"};
  for (int i=binnames.size()-1;i>=0;i--) labels.push_back(binnames[i]);


  RenameBinLabelsX(hmcinc,labels);
  RenameBinLabelsX(hmcsig,labels);
  RenameBinLabelsX(hdtinc,labels);
  RenameBinLabelsX(hmcbjt,labels);
  RenameBinLabelsX(hdtbjt,labels);
  RenameBinLabelsX(hdtb12,labels);
  RenameBinLabelsX(hmcb12,labels);
  RenameBinLabelsX(hmcb12Signal,labels);
  RenameBinLabelsX(hmcbSB,labels);

  for(unsigned i=0;i<bins.size();i++) {
    const char * end = i<bins.size()-1 ? Form("%d%d",(int)bins[i],(int)bins[i+1]) : "-1-1"; //pp==-1-1
    cout<<nbins-i+1<<endl;
    hmcinc->SetBinContent(nbins-i,res[Form("xj_mc_inc_mean%s",end)]);
    hdtinc->SetBinContent(nbins-i,res[Form("xj_data_inc_mean%s",end)]);
    hmcbjt->SetBinContent(nbins-i,res[Form("xj_mc_bjt_mean%s",end)]);
    hdtbjt->SetBinContent(nbins-i,res[Form("xj_data_bjt_mean%s",end)]);
    hdtb12->SetBinContent(nbins-i,res[Form("xj_data_b12_mean%s",end)]);
    hmcb12->SetBinContent(nbins-i,res[Form("xj_mc_b12_mean%s",end)]);
    hmcsig->SetBinContent(nbins-i,res[Form("xj_mc_sig2_inc_mean%s",end)]);
    hmcbSB->SetBinContent(nbins-i,res[Form("xj_mc_bjtSB_mean%s",end)]);
    hmcb12Signal->SetBinContent(nbins-i,res[Form("xj_mc_b12Signal_mean%s",end)]);

    hmcinc->SetBinError(nbins-i,res[Form("xj_mc_inc_meanerror%s",end)]);
    hdtinc->SetBinError(nbins-i,res[Form("xj_data_inc_meanerror%s",end)]);
    hmcbjt->SetBinError(nbins-i,res[Form("xj_mc_bjt_meanerror%s",end)]);
    hdtbjt->SetBinError(nbins-i,res[Form("xj_data_bjt_meanerror%s",end)]);
    hdtb12->SetBinError(nbins-i,res[Form("xj_data_b12_meanerror%s",end)]);
    hmcb12->SetBinError(nbins-i,res[Form("xj_mc_b12_meanerror%s",end)]);
    hmcsig->SetBinError(nbins-i,res[Form("xj_mc_sig2_inc_meanerror%s",end)]);
    hmcbSB->SetBinError(nbins-i,res[Form("xj_mc_bjtSB_meanerror%s",end)]);
    hmcb12Signal->SetBinError(nbins-i,res[Form("xj_mc_b12Signal_meanerror%s",end)]);


  }

  SetB({hmcbjt,hdtbjt,hdtb12,hmcb12,hmcb12Signal,hmcbSB});
  SetInc({hmcinc,hdtinc,hmcsig});
  SetMC({hmcinc,hmcbjt,hmcsig,hmcbSB,hmcb12,hmcb12Signal});
  SetData({hdtinc,hdtbjt,hdtb12});
  hmcsig->SetMarkerColor(kblue);
  hmcsig->SetLineColor(kblue);
  hmcbSB->SetMarkerColor(kred);
  hmcbSB->SetLineColor(kred);
  hmcsig->SetMarkerStyle(25);
  hmcbSB->SetMarkerStyle(25);
  hmcb12Signal->SetMarkerStyle(25);
  // hdtb12->SetMarkerStyle(kFullTriangleUp);
  hdtb12->SetMarkerColor(kredLight);
  hdtb12->SetLineColor(kredLight);
  hmcb12Signal->SetMarkerColor(kredLight-7);
  hmcb12Signal->SetLineColor(kredLight-7);
  hmcb12->SetMarkerColor(kredLight);
  hmcb12->SetLineColor(kredLight);

  /*
  hdtb12->GetXaxis()->SetLimits(0.1,nbins+0.1);
  hmcb12->GetXaxis()->SetLimits(0.1,nbins+0.1);
  hmcb12Signal->GetXaxis()->SetLimits(0.1,nbins+0.1);

  hmcsig->GetXaxis()->SetLimits(-0.1,nbins-0.1);
  hmcinc->GetXaxis()->SetLimits(-0.1,nbins-0.1);
  hdtinc->GetXaxis()->SetLimits(-0.1,nbins-0.1);
  */
  aktstring = "";
  plotoverwritecolors = false;
  plotlegendpos = BottomLeft;//BottomRight;
  plotymin = 0.45;
  plotymax = 0.75;

  // fitPaleLine({hmcb12Signal,hdtb12,hmcsig,hdtinc});


  hmcincsys = (TH1F *)hmcinc->Clone("hmcincsys");
  hdtincsys = (TH1F *)hdtinc->Clone("hdtincsys");
  hmcbjtsys = (TH1F *)hmcbjt->Clone("hmcbjtsys");
  hdtbjtsys = (TH1F *)hdtbjt->Clone("hdtbjtsys");

  //loadsyst();

  for(unsigned i=0;i<bins.size();i++) {
    float dx = abs(hmcinc->GetBinContent(i+1)-hmcsig->GetBinContent(i+1));
    float sig = addIn2(addIn2(hmcinc->GetBinError(i+1),hmcsig->GetBinError(i+1)),dx);
    hdtincsys->SetBinError(i+1,sig);

    dx = abs(hmcbjt->GetBinContent(i+1)-hmcbSB->GetBinContent(i+1));
    sig = addIn2(addIn2(hmcbjt->GetBinError(i+1),hmcbSB->GetBinError(i+1)),dx);
    hdtbjtsys->SetBinError(i+1,sig);
  }


hdtincsys->SetBinError(1,0.022);
hdtincsys->SetBinError(2,0.023);
hdtincsys->SetBinError(3,0.026);
hdtincsys->SetBinError(4,0.036);


hdtbjtsys->SetBinError(1,0.023);//0.023);
hdtbjtsys->SetBinError(2,0.025);//0.023);
hdtbjtsys->SetBinError(3,0.028);//0.024);
hdtbjtsys->SetBinError(4,0.035);//0.027);



hmcincsys->SetFillColor(kblueLight);   hmcincsys->SetFillStyle(1001);
hdtincsys->SetFillColor(kblueLight);   hdtincsys->SetFillStyle(1001);
hmcbjtsys->SetFillColor(kredLight);   hmcbjtsys->SetFillStyle(1001); //kRed-10); 
hdtbjtsys->SetFillColor(kredLight);   hdtbjtsys->SetFillStyle(1001); //kRed-10); 

hmcincsys->SetMarkerSize(0);
hdtincsys->SetMarkerSize(0);
hmcbjtsys->SetMarkerSize(0);
hdtbjtsys->SetMarkerSize(0);

plotputgrid = true;
Draw({hmcinc,hdtinc,hmcbjt,hdtbjt,hdtb12,hmcsig,hmcbSB,hmcb12,hmcb12Signal},"E1");

Draw({hmcinc,hdtinc,hdtb12,hmcsig,hmcb12,hmcb12Signal},"E1");

Draw({hdtb12,hmcb12,hmcb12Signal},"E1");
Draw({hmcinc,hdtinc,hmcsig},"E1");

  
TCanvas *c = new TCanvas("cmoneyinc","cmoneyinc",600,600);
TExec *ex2 = new TExec("ex2","gStyle->SetErrorX(0)");//crazy solution, bravo root
ex2->Draw();
hdtincsys->SetMinimum(0.45);
hdtincsys->SetMaximum(0.75);
hdtincsys->Draw("E2");

TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
setex1->Draw();  

hmcsig->Draw("E1,same");
hdtinc->Draw("E1,same");
TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.1)");
setex2->Draw();

   hdtinc->SetMarkerColor(kblue);
  hdtinc->SetLineColor(kblue);
 
TLegend *l = getLegend();
l->AddEntry(hdtinc,"Data","P");
l->AddEntry(hmcsig,"Simulation","P");
l->SetHeader("Inclusive dijets");
l->Draw();


 // writing the lumi information and the CMS "logo"
  CMS_lumi(c, iPeriod, iPos );

 SavePlots(c,"moneyplotINC");




TCanvas *c2 = new TCanvas("cmoneybjt","cmoneybjt",600,600);
ex2->Draw();
hdtbjtsys->SetMinimum(0.45);
hdtbjtsys->SetMaximum(0.75);
hdtbjtsys->Draw("E2");
setex1->Draw();  

hmcbSB->Draw("E1,same");
hdtbjt->Draw("E1,same");
setex2->Draw();

TLegend *l2 = getLegend();
l2->AddEntry(hdtbjt,"Data","P");
l2->AddEntry(hmcbSB,"Simulation","P");
l2->SetHeader("b-dijets");
l2->Draw();
SavePlots(c2,"moneyplotBJT");


 
TCanvas *c3 = new TCanvas("c3","c3",600,600);
 
// Transform into ncoll weighted 



 TH1F *hframe = new TH1F("hframe","hframe",1,-24.99,400);
 hframe->SetXTitle("N_{coll}-weighted N_{part}");
 hframe->SetYTitle("#LTx_{J}#GT");
 hframe->SetMaximum(0.75);
 hframe->SetMinimum(0.45);

hframe -> GetXaxis() -> CenterTitle();
hframe -> GetYaxis() -> CenterTitle();
 
 TGraphErrors *gData_PbPb_inc = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_inc_sys = new TGraphErrors(3);
 TGraphErrors *gData_pp_inc = new TGraphErrors(1);
 TGraphErrors *gData_pp_inc_sys = new TGraphErrors(1);

 TGraphErrors *gMC_PbPb_inc = new TGraphErrors(3);
 TGraphErrors *gMC_pp_inc = new TGraphErrors(1);
 
 // ncoll weighted values from Shengquan
 float npart[3] = {46.81,226.7,358.8};
 
 
 for(int i=0;i<3;i++){
   gData_PbPb_inc->SetPoint(i,npart[i],hdtinc->GetBinContent(i+2));
   gData_PbPb_inc->SetPointError(i,0.,hdtinc->GetBinError(i+2));
   
   gData_PbPb_inc_sys->SetPoint(i,npart[i],hdtinc->GetBinContent(i+2));
   gData_PbPb_inc_sys->SetPointError(i,10.,hdtincsys->GetBinError(i+2));

   gMC_PbPb_inc->SetPoint(i,npart[i],hmcsig->GetBinContent(i+2));
   gMC_PbPb_inc->SetPointError(i,0.,hmcsig->GetBinError(i+2));
 }

   gData_pp_inc->SetPoint(0,2,hdtinc->GetBinContent(1));
   gData_pp_inc->SetPointError(0,0.,hdtinc->GetBinError(1));

   gData_pp_inc_sys->SetPoint(0,2,hdtinc->GetBinContent(1));
   gData_pp_inc_sys->SetPointError(0,10.,hdtincsys->GetBinError(1));

   gMC_pp_inc->SetPoint(0,2,hmcsig->GetBinContent(1));
   gMC_pp_inc->SetPointError(0,0.,hmcsig->GetBinError(1));

 

 gData_PbPb_inc->SetMarkerColor(kblue);
 gData_PbPb_inc->SetLineColor(kblue);
 gData_PbPb_inc_sys->SetFillColor(kblueLight);

 gData_pp_inc->SetMarkerColor(kblue);
 gData_pp_inc->SetLineColor(kblue);
  gData_pp_inc_sys->SetFillColor(kblueLight);

  gMC_PbPb_inc->SetMarkerStyle(25);
  gMC_PbPb_inc->SetMarkerColor(kblue);
  gMC_PbPb_inc->SetLineColor(kblue);
  
  gMC_pp_inc->SetMarkerStyle(25);
  gMC_pp_inc->SetMarkerColor(kblue);
  gMC_pp_inc->SetLineColor(kblue);
  
 hframe->Draw();
 gData_PbPb_inc_sys->Draw("e2");
 gData_PbPb_inc->Draw("p");
 gData_pp_inc_sys->Draw("e2");
 gData_pp_inc->Draw("p");
 
 gMC_PbPb_inc->Draw("p");
 gMC_pp_inc->Draw("p");

 TText *t1a = new TText(-7,0.73,"pp");
 TText *t2a = new TText(25,0.697,"30-100%");
 TText *t3a = new TText(205,0.675,"10-30%");
 TText *t4a = new TText(340,0.665,"0-10%");

 t1a->SetTextSize(17);
 t2a->SetTextSize(17);
 t3a->SetTextSize(17);
 t4a->SetTextSize(17);
 
 t1a->Draw();
 t2a->Draw();
 t3a->Draw();
 t4a->Draw();

 
  CMS_lumi(c3, iPeriod, iPos ); 

  l->Draw();

TCanvas *c4 = new TCanvas("c4","c4",600,600);

 
 TGraphErrors *gData_PbPb_b = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_b_sys = new TGraphErrors(3);
 TGraphErrors *gData_pp_b = new TGraphErrors(1);
 TGraphErrors *gData_pp_b_sys = new TGraphErrors(1);

 TGraphErrors *gMC_PbPb_b = new TGraphErrors(3);
 TGraphErrors *gMC_pp_b = new TGraphErrors(1);
 
 
 
 for(int i=0;i<3;i++){
   gData_PbPb_b->SetPoint(i,npart[i],hdtbjt->GetBinContent(i+2));
   gData_PbPb_b->SetPointError(i,0.,hdtbjt->GetBinError(i+2));
   
   gData_PbPb_b_sys->SetPoint(i,npart[i],hdtbjt->GetBinContent(i+2));
   gData_PbPb_b_sys->SetPointError(i,10.,hdtbjtsys->GetBinError(i+2));

   gMC_PbPb_b->SetPoint(i,npart[i],hmcbSB->GetBinContent(i+2));
   gMC_PbPb_b->SetPointError(i,0.,hmcbSB->GetBinError(i+2));
 }

   gData_pp_b->SetPoint(0,2,hdtbjt->GetBinContent(1));
   gData_pp_b->SetPointError(0,0.,hdtbjt->GetBinError(1));

   gData_pp_b_sys->SetPoint(0,2,hdtbjt->GetBinContent(1));
   gData_pp_b_sys->SetPointError(0,10.,hdtbjtsys->GetBinError(1));

   gMC_pp_b->SetPoint(0,2,hmcbSB->GetBinContent(1));
   gMC_pp_b->SetPointError(0,0.,hmcbSB->GetBinError(1));

 

 gData_PbPb_b->SetMarkerColor(kred);
 gData_PbPb_b->SetLineColor(kred);
 gData_PbPb_b_sys->SetFillColor(kredLight);

 gData_pp_b->SetMarkerColor(kred);
 gData_pp_b->SetLineColor(kred);
  gData_pp_b_sys->SetFillColor(kredLight);

  gMC_PbPb_b->SetMarkerStyle(25);
  gMC_PbPb_b->SetMarkerColor(kred);
  gMC_PbPb_b->SetLineColor(kred);
  
  gMC_pp_b->SetMarkerStyle(25);
  gMC_pp_b->SetMarkerColor(kred);
  gMC_pp_b->SetLineColor(kred);
  
 hframe->Draw();
 gData_PbPb_b_sys->Draw("e2");
 gData_PbPb_b->Draw("p");
 gData_pp_b_sys->Draw("e2");
 gData_pp_b->Draw("p");
 
 gMC_PbPb_b->Draw("p");
 gMC_pp_b->Draw("p");

 TText *t1b = new TText(-7,0.702,"pp");
 TText *t2b = new TText(25,0.687,"30-100%");
 TText *t3b = new TText(205,0.665,"10-30%");
 TText *t4b = new TText(340,0.65,"0-10%");

 t1b->SetTextSize(17);
 t2b->SetTextSize(17);
 t3b->SetTextSize(17);
 t4b->SetTextSize(17);
 
 t1b->Draw();
 t2b->Draw();
 t3b->Draw();
 t4b->Draw();

 
  CMS_lumi(c4, iPeriod, iPos ); 

  l2->Draw();


  
}



  
