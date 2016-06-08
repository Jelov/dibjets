#include "../helpers/looptuple.h"
#include "../helpers/plotting.h"
#include "../helpers/physics.h"
#include "TExec.h"

TH1F *hmcincsys, *hdtincsys, *hmcbjtsys, *hdtbjtsys;

float addIn2(float v1, float v2)
{
  return sqrt(v1*v1+v2*v2);
}

// void loadsyst()
// {
//   auto mqcd = ReadFromFile("../mistag/hydjetclosureqcd/hydjetclosureqcdsyst.root");

//   hmcincsys->SetBinError(1,mqcd["closure020"]);
//   hmcincsys->SetBinError(2,mqcd["closure2060"]);
//   hmcincsys->SetBinError(3,mqcd["closure60200"]);
//   hmcincsys->SetBinError(4,0);

//   hdtincsys->SetBinError(1,mqcd["closure020"]);
//   hdtincsys->SetBinError(2,mqcd["closure2060"]);
//   hdtincsys->SetBinError(3,mqcd["closure60200"]);
//   hdtincsys->SetBinError(4,0);


//   auto sysbjt = ReadFromFile("results_0510_pthat50/results.root");
//   float bjtsyspp = addIn2(sysbjt["xjmcbjtsyslo-1-1"],sysbjt["xjmcbjtsyshi-1-1"]);
//   float bjtsys020 = addIn2(sysbjt["xjmcbjtsyslo020"],sysbjt["xjmcbjtsyshi020"]);
//   float bjtsys2060 = addIn2(sysbjt["xjmcbjtsyslo2060"],sysbjt["xjmcbjtsyshi2060"]);
//   float bjtsys60200 = addIn2(sysbjt["xjmcbjtsyslo60200"],sysbjt["xjmcbjtsyshi60200"]);


//   auto mbjt = ReadFromFile("../mistag/hydjetclosure_upd/hydjetclosuresyst.root");

//   float bjtfullsyspp = bjtsyspp;
//   float bjtfullsys020 = addIn2(bjtsys020,mbjt["closure020"]);
//   float bjtfullsys2060 = addIn2(bjtsys2060,mbjt["closure2060"]);
//   float bjtfullsys60200 = addIn2(bjtsys60200,mbjt["closure60200"]);


//   hmcbjtsys->SetBinError(1,bjtfullsys020);
//   hmcbjtsys->SetBinError(2,bjtfullsys2060);
//   hmcbjtsys->SetBinError(3,bjtfullsys60200);
//   hmcbjtsys->SetBinError(4,bjtfullsyspp);

//   hdtbjtsys->SetBinError(1,bjtfullsys020);
//   hdtbjtsys->SetBinError(2,bjtfullsys2060);
//   hdtbjtsys->SetBinError(3,bjtfullsys60200);
//   hdtbjtsys->SetBinError(4,bjtfullsyspp);



// }

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
  macro m("moneyplot_"+name);

  auto res = ReadFromFile("results_"+name+"/results.root");
  for (auto i:res) {
    cout<<i.first<<" = "<<i.second<<endl;
  }

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
    hmcinc->SetBinContent(nbins-i,res[Form("xj_mc_inc_mean%s",end)]);      //i+1
    hdtinc->SetBinContent(nbins-i,res[Form("xj_data_inc_mean%s",end)]);      //i+1
    hmcbjt->SetBinContent(nbins-i,res[Form("xj_mc_bjt_mean%s",end)]);      //i+1
    hdtbjt->SetBinContent(nbins-i,res[Form("xj_data_bjt_mean%s",end)]);      //i+1
    hdtb12->SetBinContent(nbins-i,res[Form("xj_data_b12_mean%s",end)]);      //i+1
    hmcb12->SetBinContent(nbins-i,res[Form("xj_mc_b12_mean%s",end)]);      //i+1
    hmcsig->SetBinContent(nbins-i,res[Form("xj_mc_sig2_inc_mean%s",end)]);      //i+1
    hmcbSB->SetBinContent(nbins-i,res[Form("xj_mc_bjtSB_mean%s",end)]);      //i+1
    hmcb12Signal->SetBinContent(nbins-i,res[Form("xj_mc_b12Signal_mean%s",end)]);      //i+1

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
  hmcsig->SetMarkerColor(kBlue-7);
  hmcsig->SetLineColor(kBlue-7);
  hmcbSB->SetMarkerColor(kRed-7);
  hmcbSB->SetLineColor(kRed-7);
  hmcsig->SetMarkerStyle(kOpenDiamond);
  hmcbSB->SetMarkerStyle(kOpenDiamond);
  hmcb12Signal->SetMarkerStyle(kOpenDiamond);
  // hdtb12->SetMarkerStyle(kFullTriangleUp);
  hdtb12->SetMarkerColor(kMagenta);
  hdtb12->SetLineColor(kMagenta);
  hmcb12Signal->SetMarkerColor(kMagenta-7);
  hmcb12Signal->SetLineColor(kMagenta-7);
  hmcb12->SetMarkerColor(kMagenta);
  hmcb12->SetLineColor(kMagenta);

  hdtb12->GetXaxis()->SetLimits(0.1,nbins+0.1);
  hmcb12->GetXaxis()->SetLimits(0.1,nbins+0.1);
  hmcb12Signal->GetXaxis()->SetLimits(0.1,nbins+0.1);

  hmcsig->GetXaxis()->SetLimits(-0.1,nbins-0.1);
  hmcinc->GetXaxis()->SetLimits(-0.1,nbins-0.1);
  hdtinc->GetXaxis()->SetLimits(-0.1,nbins-0.1);

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



hmcincsys->SetFillColor(kCyan-10);   hmcincsys->SetFillStyle(1001);
hdtincsys->SetFillColor(kCyan-10);   hdtincsys->SetFillStyle(1001);
hmcbjtsys->SetFillColor(kRed-10);   hmcbjtsys->SetFillStyle(1001); //kRed-10); 
hdtbjtsys->SetFillColor(kRed-10);   hdtbjtsys->SetFillStyle(1001); //kRed-10); 

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

TLegend *l = getLegend();
l->AddEntry(hmcsig,"MC","P");
l->AddEntry(hdtinc,"Data","P");
l->SetHeader("Inclusive dijets");
l->Draw();
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
l2->AddEntry(hmcbSB,"MC","P");
l2->AddEntry(hdtbjt,"Data","P");
l2->SetHeader("b-dijets");
l2->Draw();
SavePlots(c2,"moneyplotBJT");


}



  
// TCanvas *c = new TCanvas("cmoney","cmoney",600,600);
// TExec *ex2 = new TExec("ex2","gStyle->SetErrorX(0)");//crazy solution, bravo root
// ex2->Draw();

// hdtincsys->SetMinimum(0.45);//0.55
// hdtincsys->SetMaximum(0.75);
// //hmcincsys->Draw("E2");
// hdtincsys->Draw("E2");
// //hmcbjtsys->Draw("E2,same");
// hdtbjtsys->Draw("E2,same");

// TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
// setex1->Draw();  

// //hmcinc->Draw("E1,same");
// hmcsig->Draw("E1,same");
// hdtinc->Draw("E1,same");
// //hmcbjt->Draw("E1,same");
// hmcbSB->Draw("E1,same");
// hdtbjt->Draw("E1,same");
// TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.1)");
// setex2->Draw();

// TLegend *l = getLegend();
// l->AddEntry(hmcinc,hmcinc->GetTitle(),"P");
// l->AddEntry(hdtinc,hdtinc->GetTitle(),"P");
// l->AddEntry(hmcbjt,hmcbjt->GetTitle(),"P");
// l->AddEntry(hdtbjt,hdtbjt->GetTitle(),"P");
// l->Draw();
// SavePlots(c,"moneyplot");


