#include "../helpers/looptuple.h"
#include "../helpers/plotting.h"
#include "../helpers/physics.h"
// #include "./CMS_lumi.C"
#include "TExec.h"

 // ncoll weighted values from Shengquan
 float npart[3] = {102.,239.9,363.4};
 

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

void CopyToGraph(TGraphErrors *g, TH1F *h)
{
    for(int i=0;i<3;i++){
    g->SetPoint(i+1,npart[i],h->GetBinContent(i+2));
    g->SetPointError(i+1,0.,h->GetBinError(i+2));
  }
  g->SetPoint(0,2,h->GetBinContent(1));
  g->SetPointError(0,0.,h->GetBinError(1));
  cout<<"???"<<h->GetBinContent(1)<<endl;
}

void moneyplot(TString name="")
{
  bool drawcalo = false;

  bool drawsimulation = false;
  bool drawsmeareddata = true;
  bool drawsmearedsys = true;
  bool drawsmearedmc = false;
  bool drawsys = true;

  // Lumi stuff

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = lumi_sqrtSppPbPb;       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  //int iPos=11; // : top-left, left-aligned
  int  iPos=33;// : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)


  macro m("moneyplot_"+name);

  auto res = ReadFromFile("results_"+name+"_default/results.root");
  auto resppsmeared1 = ReadFromFile("results_"+name+"_ppsmear1/results.root");
  auto resppsmeared2 = ReadFromFile("results_"+name+"_ppsmear2/results.root");
  auto resppsmeared3 = ReadFromFile("results_"+name+"_ppsmear3/results.root");

  //  auto calores = ReadFromFile("results_0704_calotrig/results.root");

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

  auto hdtb12calo   = geth("hdtb12calo","Data b-jets 12 Calo;;#LTx_{J}#GT");

  auto hdtipp       = geth("hdtpps","Data pp Smeared;;#LTx_{J}#GT");
  auto hdtbpp       = geth("hdtbpp","Data b-jets Smeared;;#LTx_{J}#GT");
  auto hmcipp       = geth("hmcipp","MC qcd Smeared;;#LTx_{J}#GT");
  auto hmcbpp       = geth("hmcbpp","MC b-jet Smeared;;#LTx_{J}#GT");

  
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
  RenameBinLabelsX(hdtipp,labels);
  RenameBinLabelsX(hdtbpp,labels);
  RenameBinLabelsX(hmcipp,labels);
  RenameBinLabelsX(hmcbpp,labels);


  for(unsigned i=0;i<bins.size();i++) {
    const char * end = i<bins.size()-1 ? Form("%d%d",(int)bins[i],(int)bins[i+1]) : "pp"; //pp==-1-1
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
    //    hdtb12calo->SetBinContent(nbins-i,calores[Form("xj_data_b12_mean%s",end)]);

    hmcinc->SetBinError(nbins-i,res[Form("xj_mc_inc_meanerror%s",end)]);
    hdtinc->SetBinError(nbins-i,res[Form("xj_data_inc_meanerror%s",end)]);
    hmcbjt->SetBinError(nbins-i,res[Form("xj_mc_bjt_meanerror%s",end)]);
    hdtbjt->SetBinError(nbins-i,res[Form("xj_data_bjt_meanerror%s",end)]);
    hdtb12->SetBinError(nbins-i,res[Form("xj_data_b12_meanerror%s",end)]);
    hmcb12->SetBinError(nbins-i,res[Form("xj_mc_b12_meanerror%s",end)]);
    hmcsig->SetBinError(nbins-i,res[Form("xj_mc_sig2_inc_meanerror%s",end)]);
    hmcbSB->SetBinError(nbins-i,res[Form("xj_mc_bjtSB_meanerror%s",end)]);
    hmcb12Signal->SetBinError(nbins-i,res[Form("xj_mc_b12Signal_meanerror%s",end)]);
    //    hdtb12calo->SetBinError(nbins-i,calores[Form("xj_data_b12_meanerror%s",end)]);

  }

  hdtipp->SetBinContent(nbins-0,resppsmeared1["xj_data_inc_meanpp"]);
  hdtipp->SetBinError(nbins-0,resppsmeared1["xj_data_inc_meanerrorpp"]);
  hdtbpp->SetBinContent(nbins-0,resppsmeared1["xj_data_b12_meanpp"]);
  hdtbpp->SetBinError(nbins-0,resppsmeared1["xj_data_b12_meanerrorpp"]);
  hmcipp->SetBinContent(nbins-0,resppsmeared1["xj_mc_inc_meanpp"]);
  hmcipp->SetBinError(nbins-0,resppsmeared1["xj_mc_inc_meanerrorpp"]);
  hmcbpp->SetBinContent(nbins-0,resppsmeared1["xj_mc_b12Signal_meanpp"]);
  hmcbpp->SetBinError(nbins-0,resppsmeared1["xj_mc_b12Signal_meanerrorpp"]);

  hdtipp->SetBinContent(nbins-1,resppsmeared2["xj_data_inc_meanpp"]);
  hdtipp->SetBinError(nbins-1,resppsmeared2["xj_data_inc_meanerrorpp"]);
  hdtbpp->SetBinContent(nbins-1,resppsmeared2["xj_data_b12_meanpp"]);
  hdtbpp->SetBinError(nbins-1,resppsmeared2["xj_data_b12_meanerrorpp"]);
  hmcipp->SetBinContent(nbins-1,resppsmeared2["xj_mc_inc_meanpp"]);
  hmcipp->SetBinError(nbins-1,resppsmeared2["xj_mc_inc_meanerrorpp"]);
  hmcbpp->SetBinContent(nbins-1,resppsmeared2["xj_mc_b12Signal_meanpp"]);
  hmcbpp->SetBinError(nbins-1,resppsmeared2["xj_mc_b12Signal_meanerrorpp"]);

  hdtipp->SetBinContent(nbins-2,resppsmeared3["xj_data_inc_meanpp"]);
  hdtipp->SetBinError(nbins-2,resppsmeared3["xj_data_inc_meanerrorpp"]);
  hdtbpp->SetBinContent(nbins-2,resppsmeared3["xj_data_b12_meanpp"]);
  hdtbpp->SetBinError(nbins-2,resppsmeared3["xj_data_b12_meanerrorpp"]);
  hmcipp->SetBinContent(nbins-2,resppsmeared3["xj_mc_inc_meanpp"]);
  hmcipp->SetBinError(nbins-2,resppsmeared3["xj_mc_inc_meanerrorpp"]);
  hmcbpp->SetBinContent(nbins-2,resppsmeared3["xj_mc_b12Signal_meanpp"]);
  hmcbpp->SetBinError(nbins-2,resppsmeared3["xj_mc_b12Signal_meanerrorpp"]);

  hmcipp->SetMarkerStyle(kOpenCircle);
  hmcbpp->SetMarkerStyle(kOpenCircle);

  SetB({hmcbjt,hdtbjt,hdtb12,hmcb12,hmcb12Signal,hmcbSB,hdtbpp});//,hdtbpp
  SetInc({hmcinc,hdtinc,hdtipp,hmcsig});
  SetMC({hmcinc,hmcbjt,hmcsig,hmcbSB,hmcb12,hmcb12Signal,hdtipp,hdtbpp});
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

  //  hdtb12calo->SetMarkerStyle(24);
  //  hdtb12calo->SetMarkerColor(kRed+3);
  //  hdtb12calo->SetLineColor(kRed+3);


  hdtipp->SetMarkerColor(kblue);

// hdtbpp

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
  plotymin = 0.55;
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


hdtbjtsys->SetBinError(1,0.008);
hdtbjtsys->SetBinError(2,0.014);
hdtbjtsys->SetBinError(3,0.018);
hdtbjtsys->SetBinError(4,0.023);


hdtincsys->SetBinError(1,0.007);
hdtincsys->SetBinError(2,0.010);
hdtincsys->SetBinError(3,0.016);
hdtincsys->SetBinError(4,0.023);


hmcincsys->SetFillColor(kblueLight);   hmcincsys->SetFillStyle(1001);
hdtincsys->SetFillColor(kblueLight);   hdtincsys->SetFillStyle(1001);
hmcbjtsys->SetFillColor(kredLight);   hmcbjtsys->SetFillStyle(1001); //kRed-10); 
hdtbjtsys->SetFillColor(kredLight);   hdtbjtsys->SetFillStyle(1001); //kRed-10); 

hmcincsys->SetMarkerSize(0);
hdtincsys->SetMarkerSize(0);
hmcbjtsys->SetMarkerSize(0);
hdtbjtsys->SetMarkerSize(0);

plotputgrid = true;
// Draw({hmcinc,hdtinc,hmcbjt,hdtbjt,hdtb12,hmcsig,hmcbSB,hmcb12,hmcb12Signal},"E1");

// Draw({hmcinc,hdtinc,hdtb12,hmcsig,hmcb12,hmcb12Signal},"E1");

// Draw({hdtb12,hmcb12,hmcb12Signal},"E1");
// Draw({hmcinc,hdtinc,hmcsig},"E1");

  
TCanvas *c = new TCanvas("cmoneyinc","cmoneyinc",600,600);
TExec *ex2 = new TExec("ex2","gStyle->SetErrorX(0)");//crazy solution, bravo root
ex2->Draw();
hdtincsys->SetMinimum(0.55);
hdtincsys->SetMaximum(0.75);
hdtincsys->Draw("E2");

TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
setex1->Draw();  

hmcsig->Draw("E1,same");
hdtinc->Draw("E1,same");

hdtipp->Draw("E1,same");

TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.1)");
setex2->Draw();

   hdtinc->SetMarkerColor(kblue);
  hdtinc->SetLineColor(kblue);
 
TLegend *l = getLegend();
l->AddEntry(hdtinc,"Data","P");
if (drawsimulation) l->AddEntry(hmcsig,"Simulation","P");
if (drawsmeareddata) l->AddEntry(hdtipp,drawsmeareddata ? "pp-based reference" : "Simulation","P");
if (drawsmearedmc) l->AddEntry(hmcipp,"Pythia smeared","P");
l->SetHeader("Inclusive dijets");
l->Draw();


 // writing the lumi information and the CMS "logo"
  CMS_lumi(c, iPeriod, iPos );

 // SavePlot(c,"moneyplotINC");




TCanvas *c2 = new TCanvas("cmoneybjt","cmoneybjt",600,600);
ex2->Draw();
hdtbjtsys->SetMinimum(0.55);
hdtbjtsys->SetMaximum(0.75);
hdtbjtsys->Draw("E2");
setex1->Draw();  

hmcbSB->Draw("E1,same");
hdtbjt->Draw("E1,same");

hdtbpp->Draw("E1,same");

setex2->Draw();

TLegend *l2 = getLegend();
l2->AddEntry(hdtbjt,"Data","P");
//if (drawcalo) l2->AddEntry(hdtb12calo,"Data CaloJet trigger","P");
if (drawsimulation) l2->AddEntry(hmcbSB,"Simulation","P");
if (drawsmeareddata) l2->AddEntry(hdtbpp,drawsmeareddata ? "pp-based reference" : "Simulation","P");
if (drawsmearedmc) l2->AddEntry(hmcbpp,"Pythia smeared","P");
l2->SetHeader("b dijets");
l2->Draw();
// SavePlot(c2,"moneyplotBJT");


 
TCanvas *c3 = new TCanvas("c3","c3",600,600);
 
// Transform into ncoll weighted 



 TH1F *hframe = new TH1F("hframe","hframe",1,-24.99,400);
 hframe->SetXTitle("#LTN_{part}#GT (N_{coll}-weighted)");
 hframe->SetYTitle("#LTx_{J}#GT");
 hframe->SetMaximum(0.75);
 hframe->SetMinimum(0.55);

hframe -> GetXaxis() -> CenterTitle();
hframe -> GetYaxis() -> CenterTitle();
 
 TGraphErrors *gData_PbPb_inc = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_inc_sys = new TGraphErrors(3);
 TGraphErrors *gData_pp_inc = new TGraphErrors(1);
 TGraphErrors *gData_pp_inc_sys = new TGraphErrors(1);

 TGraphErrors *gMC_PbPb_inc = new TGraphErrors(3);
 TGraphErrors *gMC_pp_inc = new TGraphErrors(1);

 TGraphErrors *gData_PbPb_ppsmeared = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_qcdsmeared = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_bjtsmeared = new TGraphErrors(3);

 TGraphErrors *gData_PbPb_ppsmearedsys = new TGraphErrors(3);

 

 
 for(int i=0;i<3;i++){
   gData_PbPb_inc->SetPoint(i,npart[i],hdtinc->GetBinContent(i+2));
   gData_PbPb_inc->SetPointError(i,0.,hdtinc->GetBinError(i+2));
   
   gData_PbPb_inc_sys->SetPoint(i,npart[i],hdtinc->GetBinContent(i+2));
   gData_PbPb_inc_sys->SetPointError(i,10.,hdtincsys->GetBinError(i+2));

   gMC_PbPb_inc->SetPoint(i,npart[i],hmcsig->GetBinContent(i+2));
   gMC_PbPb_inc->SetPointError(i,0.,hmcsig->GetBinError(i+2));

   gData_PbPb_ppsmeared->SetPoint(i,npart[i],hdtipp->GetBinContent(i+2));
   gData_PbPb_ppsmeared->SetPointError(i,0.,hdtipp->GetBinError(i+2));

   gData_PbPb_ppsmearedsys->SetPoint(i,npart[i],hdtipp->GetBinContent(i+2));
   gData_PbPb_ppsmearedsys->SetPointError(i,10.,hdtincsys->GetBinError(1));

   gData_PbPb_qcdsmeared->SetPoint(i,npart[i],hmcipp->GetBinContent(i+2)); 
   gData_PbPb_qcdsmeared->SetPointError(i,0.,hmcipp->GetBinError(i+2));

   gData_PbPb_bjtsmeared->SetPoint(i,npart[i],hmcbpp->GetBinContent(i+2)); 
   gData_PbPb_bjtsmeared->SetPointError(i,0.,hmcbpp->GetBinError(i+2));
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
  
// gData_PbPb_qcdsmeared->SetMarkerStyle(kOpenCircle);
// gData_PbPb_bjtsmeared->SetMarkerStyle(kOpenCircle);

  gData_PbPb_ppsmearedsys->SetFillStyle(3005);
  gData_PbPb_ppsmearedsys->SetFillColor(kblue);

    gData_PbPb_ppsmeared->SetMarkerColor(kblue);
    gData_PbPb_ppsmeared->SetLineColor(kblue);
  gData_PbPb_ppsmeared->SetMarkerStyle(kOpenSquare);

 hframe->Draw();


 if (drawsys) gData_PbPb_inc_sys->Draw("e2");
 gData_PbPb_inc->Draw("p");
 if (drawsys) gData_pp_inc_sys->Draw("e2");
 gData_pp_inc->Draw("p");

if (drawsmeareddata) gData_PbPb_ppsmeared->Draw("p");
if (drawsmearedmc) gData_PbPb_qcdsmeared->Draw("p");

if (drawsmearedsys) gData_PbPb_ppsmearedsys->Draw("e2");
 
 if (drawsimulation) gMC_PbPb_inc->Draw("p");
 if (drawsimulation) gMC_pp_inc->Draw("p");



 TText *t1a = new TText(-7,0.72,"pp"); //-7,0.73,
 TText *t2a = new TText(75,0.705,"30-100%"); //100,0.697,
 TText *t3a = new TText(215,0.69,"10-30%"); //205,0.675,
 TText *t4a = new TText(340,0.680,"0-10%"); //340,0.665,

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

  SavePlot(c3,Form("moneyincnice_%s%s",drawsmeareddata ? "pp" : "mc", drawsmearedsys ? "sys" : ""));

TCanvas *c4 = new TCanvas("c4","c4",600,600);

 
 TGraphErrors *gData_PbPb_b = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_b_sys = new TGraphErrors(3);
 TGraphErrors *gData_pp_b = new TGraphErrors(1);
 TGraphErrors *gData_pp_b_sys = new TGraphErrors(1);

 TGraphErrors *gMC_PbPb_b = new TGraphErrors(3);
 TGraphErrors *gMC_pp_b = new TGraphErrors(1);

 TGraphErrors *gData_PbPb_b_calo = new TGraphErrors(3);
 
 TGraphErrors *gData_PbPb_b_ppsmeared = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_b_ppsmearedsys = new TGraphErrors(3);
 
 
 for(int i=0;i<3;i++){
   gData_PbPb_b->SetPoint(i,npart[i],hdtb12->GetBinContent(i+2));
   gData_PbPb_b->SetPointError(i,0.,hdtb12->GetBinError(i+2));
   
   gData_PbPb_b_sys->SetPoint(i,npart[i],hdtb12->GetBinContent(i+2));
   gData_PbPb_b_sys->SetPointError(i,10.,hdtbjtsys->GetBinError(i+2));

   gMC_PbPb_b->SetPoint(i,npart[i],hmcb12Signal->GetBinContent(i+2));
   gMC_PbPb_b->SetPointError(i,0.,hmcb12Signal->GetBinError(i+2));

   //   gData_PbPb_b_calo->SetPoint(i,npart[i],hdtb12calo->GetBinContent(i+2));
   //   gData_PbPb_b_calo->SetPointError(i,0.,hdtb12calo->GetBinError(i+2));

   gData_PbPb_b_ppsmeared->SetPoint(i,npart[i],hdtbpp->GetBinContent(i+2));
   gData_PbPb_b_ppsmeared->SetPointError(i,0.,hdtbpp->GetBinError(i+2));

   gData_PbPb_b_ppsmearedsys->SetPoint(i,npart[i],hdtbpp->GetBinContent(i+2));
   gData_PbPb_b_ppsmearedsys->SetPointError(i,10.,hdtbjtsys->GetBinError(1));

 }

   gData_pp_b->SetPoint(0,2,hdtb12->GetBinContent(1));
   gData_pp_b->SetPointError(0,0.,hdtb12->GetBinError(1));

   gData_pp_b_sys->SetPoint(0,2,hdtb12->GetBinContent(1));
   gData_pp_b_sys->SetPointError(0,10.,hdtbjtsys->GetBinError(1));

   gMC_pp_b->SetPoint(0,2,hmcb12Signal->GetBinContent(1));
   gMC_pp_b->SetPointError(0,0.,hmcb12Signal->GetBinError(1));

 

 gData_PbPb_b->SetMarkerColor(kred);
 gData_PbPb_b->SetLineColor(kred);
 gData_PbPb_b_sys->SetFillColor(kredLight);

  gData_PbPb_b_calo->SetMarkerStyle(24);
 gData_PbPb_b_calo->SetMarkerColor(kRed+3);
 gData_PbPb_b_calo->SetLineColor(kRed+3);

 gData_pp_b->SetMarkerColor(kred);
 gData_pp_b->SetLineColor(kred);
  gData_pp_b_sys->SetFillColor(kredLight);

  gMC_PbPb_b->SetMarkerStyle(25);
  gMC_PbPb_b->SetMarkerColor(kred);
  gMC_PbPb_b->SetLineColor(kred);
  
  gMC_pp_b->SetMarkerStyle(25);
  gMC_pp_b->SetMarkerColor(kred);
  gMC_pp_b->SetLineColor(kred);
  
  gData_PbPb_b_ppsmearedsys->SetFillStyle(3005);
  gData_PbPb_b_ppsmearedsys->SetFillColor(kred);

    gData_PbPb_b_ppsmeared->SetMarkerColor(kred);
    gData_PbPb_b_ppsmeared->SetLineColor(kred);
  gData_PbPb_b_ppsmeared->SetMarkerStyle(kOpenSquare);
 hframe->Draw();

 if (drawsys)  gData_PbPb_b_sys->Draw("e2");
 gData_PbPb_b->Draw("p");
 if (drawsys)  gData_pp_b_sys->Draw("e2");
 gData_pp_b->Draw("p");
 if (drawcalo) gData_PbPb_b_calo->Draw("p");

 if (drawsmeareddata) gData_PbPb_b_ppsmeared->Draw("p");
 if (drawsmearedmc) gData_PbPb_bjtsmeared->Draw("p");
 
 if (drawsmearedsys) gData_PbPb_b_ppsmearedsys->Draw("e2");

 if (drawsimulation) gMC_PbPb_b->Draw("p");
 if (drawsimulation) gMC_pp_b->Draw("p");



 TText *t1b = new TText(-7,0.72,"pp"); //-7,0.702
 TText *t2b = new TText(75,0.705,"30-100%"); //25,0.687
 TText *t3b = new TText(215,0.69,"10-30%"); //205,0.665
 TText *t4b = new TText(340,0.680,"0-10%"); //340,0.65

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

  SavePlot(c4,Form("moneybjtnice_%s%s",drawsmeareddata ? "pp" : "mc", drawsmearedsys ? "sys" : ""));

  auto hbjtincratio = (TH1F *)hdtb12->Clone("hbjtincratio");
  auto hbjtincratiomc = (TH1F *)hmcb12->Clone("hbjtincratio");
  auto hincdtmcratio = (TH1F *)hdtb12->Clone("hincdtmcratio");
    
  auto hbjtincratiosys = (TH1F *)hdtbjtsys->Clone("hbjtincratiosys");
   
  hbjtincratio->Divide(hdtb12,hdtinc);
  hbjtincratiomc->Divide(hmcb12Signal,hmcsig);
  hincdtmcratio->Divide(hmcsig,hdtinc);    


  auto hbjtPbppratio = (TH1F *)hdtb12->Clone("hbjtPbppratio");
  auto hincPbppratio = (TH1F *)hdtb12->Clone("hincPbppratio");

  auto hbjtincppratio = (TH1F *)hdtb12->Clone("hbjtincppratio");

bool Pbppratio = false;

  if (Pbppratio) {
    hbjtPbppratio->Divide(hdtb12,hdtbpp);
    hincPbppratio->Divide(hdtinc,hdtipp);
  } else {
    hbjtPbppratio->Add(hdtb12,hdtbpp,1,-1);
    hincPbppratio->Add(hdtinc,hdtipp,1,-1);
  }


  hbjtincppratio->Divide(hdtbpp,hdtipp);








hbjtincratiosys->SetBinError(1,0.004);
hbjtincratiosys->SetBinError(2,0.015);
hbjtincratiosys->SetBinError(3,0.020);
hbjtincratiosys->SetBinError(4,0.032);

 
 TGraphErrors *gData_PbPb_r = new TGraphErrors(3);
 TGraphErrors *gData_PbPb_r_sys = new TGraphErrors(3);
 TGraphErrors *gData_pp_r = new TGraphErrors(1);
 TGraphErrors *gData_pp_r_sys = new TGraphErrors(1);
 TGraphErrors *gMC_PbPb_r = new TGraphErrors(3);
 TGraphErrors *gMC_pp_r = new TGraphErrors(1);

TGraphErrors *gMC_r = new TGraphErrors(4);
TGraphErrors *gDataMC_r = new TGraphErrors(4);
 
 
  TGraphErrors *gppSmeared_r = new TGraphErrors(4);


  for(int i=0;i<3;i++){
    gData_PbPb_r->SetPoint(i,npart[i],hbjtincratio->GetBinContent(i+2));
    gData_PbPb_r->SetPointError(i,0.,hbjtincratio->GetBinError(i+2));
    
    gData_PbPb_r_sys->SetPoint(i,npart[i],hbjtincratio->GetBinContent(i+2));
    gData_PbPb_r_sys->SetPointError(i,10.,hbjtincratiosys->GetBinError(i+2));

    gMC_r->SetPoint(i+1,npart[i],hbjtincratiomc->GetBinContent(i+2));
    gMC_r->SetPointError(i+1,0.,hbjtincratiomc->GetBinError(i+2));

    gDataMC_r->SetPoint(i+1,npart[i],hincdtmcratio->GetBinContent(i+2)); 
    gDataMC_r->SetPointError(i+1,0.,hincdtmcratio->GetBinError(i+2));

    gppSmeared_r->SetPoint(i+1,npart[i],hbjtincppratio->GetBinContent(i+2)); 
    gppSmeared_r->SetPointError(i+1,0.,hbjtincppratio->GetBinError(i+2));
  }

    gData_pp_r->SetPoint(0,2,hbjtincratio->GetBinContent(1));
    gData_pp_r->SetPointError(0,0.,hbjtincratio->GetBinError(1));
  
    gData_pp_r_sys->SetPoint(0,2,hbjtincratio->GetBinContent(1));
    gData_pp_r_sys->SetPointError(0,10.,hbjtincratiosys->GetBinError(1));
  
  gDataMC_r->SetPoint(0,2,hincdtmcratio->GetBinContent(1));
  gDataMC_r->SetPointError(0,0.,hincdtmcratio->GetBinError(1));

  gMC_r->SetPoint(0,2,hbjtincratiomc->GetBinContent(1));
  gMC_r->SetPointError(0,0.,hbjtincratiomc->GetBinError(1));

  gppSmeared_r->SetPoint(0,2,hbjtincppratio->GetBinContent(1));
  gppSmeared_r->SetPointError(0,0.,hbjtincppratio->GetBinError(1));


gMC_r->SetMarkerColor(kBlue);
gMC_r->SetLineColor(kBlue);
gMC_r->SetMarkerStyle(kOpenSquare);

gppSmeared_r->SetMarkerColor(kBlue);
gppSmeared_r->SetLineColor(kBlue);
gppSmeared_r->SetMarkerStyle(kOpenSquare);


gDataMC_r->SetMarkerColor(kBlack);

  gData_PbPb_r->SetMarkerColor(kRed);
  gData_PbPb_r->SetLineColor(kRed);
  gData_PbPb_r_sys->SetFillColor(kYellow);
  
  gData_pp_r->SetMarkerColor(kRed);
  gData_pp_r->SetLineColor(kRed);
  gData_pp_r_sys->SetFillColor(kYellow);

  // float xmin = gData_PbPb_r->GetXaxis()->GetMinimum();
  // float xmax = gData_PbPb_r->GetXaxis()->GetMaximum();


  hframe->SetMinimum(0.9);
  hframe->SetMaximum(1.15);
  hframe->SetYTitle("b-jet #LTx_{J}#GT / inclusive jet #LTx_{J}#GT");



  TCanvas *c5 = new TCanvas("c4","c4",600,600);

 hframe->Draw();

  TLine *line = new TLine(-7,1,400,1);
  line->SetLineColor(kGray);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw();


 gData_PbPb_r_sys->Draw("e2");
 gData_PbPb_r->Draw("p");
 gData_pp_r_sys->Draw("e2");
 gData_pp_r->Draw("p");
 
 gMC_PbPb_r->Draw("p");
 gMC_pp_r->Draw("p");
gMC_r->Draw("p");
// gppSmeared_r->Draw("p");
gDataMC_r->Draw("p");

  auto lr = new TLegend(0.2,0.2,0.4,0.35);
  lr->AddEntry(gData_PbPb_r,"Data","P");
  lr->AddEntry(gMC_r,"Pythia 6 + Hydjet","P");
  // lr->AddEntry(gppSmeared_r,"pp-based reference","P");
  
  lr->AddEntry(gDataMC_r,"Inclusive MC/Data","P");
  lr->Draw();

 TText *t1r = new TText(-7,1.015,"pp"); //-7,0.702
 TText *t2r = new TText(75,1.05,"30-100%"); //25,0.687
 TText *t3r = new TText(215,1.06,"10-30%"); //205,0.665
 TText *t4r = new TText(340,1.03,"0-10%"); //340,0.65

 t1r->SetTextSize(17);
 t2r->SetTextSize(17);
 t3r->SetTextSize(17);
 t4r->SetTextSize(17);
 
 t1r->Draw();
 t2r->Draw();
 t3r->Draw();
 t4r->Draw();

  CMS_lumi(c5, iPeriod, iPos ); 

SavePlot(c5,"bjtincratio_withunqueched");

auto hdtbppsys = (TH1F *)hdtbpp->Clone("hdtbppsys");
auto hdtippsys = (TH1F *)hdtipp->Clone("hdtippsys");

for (int i=1;i<=hdtbjtsys->GetNbinsX();i++) {
  hdtbjtsys->SetBinContent(i,hdtb12->GetBinContent(i));
  hdtincsys->SetBinContent(i,hdtinc->GetBinContent(i));

  hdtbppsys->SetBinError(i,hdtbjtsys->GetBinError(1));
  hdtippsys->SetBinError(i,hdtincsys->GetBinError(1));
}

if (Pbppratio) {
  hdtbjtsys->Divide(hdtbppsys);
  hdtincsys->Divide(hdtippsys);
} else {
  hdtbjtsys->Add(hdtbppsys,-1);
  hdtincsys->Add(hdtippsys,-1);  
}

  for (int i=0;i<3;i++) {
    gData_PbPb_b_sys->SetPoint(i+2,npart[i],hdtbjtsys->GetBinContent(i+2));
    gData_PbPb_b_sys->SetPointError(i+2,10,hdtbjtsys->GetBinError(i+2));

    gData_PbPb_inc_sys->SetPoint(i+2,npart[i],hdtincsys->GetBinContent(i+2));
    gData_PbPb_inc_sys->SetPointError(i+2,10,hdtincsys->GetBinError(i+2));

  }

gData_PbPb_b_sys->SetFillColorAlpha(gData_PbPb_b_sys->GetFillColor(),0.5);
gData_PbPb_inc_sys->SetFillColorAlpha(gData_PbPb_inc_sys->GetFillColor(),0.5);

gData_PbPb_b_sys->SetLineColor(kred);
gData_PbPb_inc_sys->SetLineColor(kblue);

gData_PbPb_b_sys->SetLineWidth(2);
gData_PbPb_inc_sys->SetLineWidth(2);





 TGraphErrors *gData_bjt_r = new TGraphErrors(4);
 TGraphErrors *gData_inc_r = new TGraphErrors(4);

 // gData_bjt_r->SetPoint(0,2,1);
 // gData_bjt_r->SetPointError(0,2,0);
 // gData_inc_r->SetPoint(0,2,1);
 // gData_inc_r->SetPointError(0,2,0);



  // for(int i=0;i<3;i++){
  //   gData_bjt_r->SetPoint(i+1,npart[i],hbjtPbppratio->GetBinContent(i+2));
  //   gData_bjt_r->SetPointError(i+1,0.,hbjtPbppratio->GetBinError(i+2));
  //   gData_inc_r->SetPoint(i+1,npart[i],hincPbppratio->GetBinContent(i+2));
  //   gData_inc_r->SetPointError(i+1,0.,hincPbppratio->GetBinError(i+2));
  // }
CopyToGraph(gData_bjt_r,hbjtPbppratio);
CopyToGraph(gData_inc_r,hincPbppratio);

  gData_inc_r->SetMarkerStyle(kFullCircle);
  gData_inc_r->SetMarkerColor(kblue);
  gData_inc_r->SetLineColor(kblue);

  gData_bjt_r->SetMarkerStyle(kFullCircle);
  gData_bjt_r->SetMarkerColor(kred);
  gData_bjt_r->SetLineColor(kred);



  if (Pbppratio) {
    hframe->SetMinimum(0.85);
    hframe->SetMaximum(1.05);
    hframe->SetYTitle("PbPb #LTx_{J}#GT / pp #LTx_{J}#GT");
  } else {
    hframe->SetMinimum(-0.15);
    hframe->SetMaximum(0.05);
    hframe->SetYTitle("PbPb #LTx_{J}#GT - pp #LTx_{J}#GT");

  }



  TCanvas *c6 = new TCanvas("c4","c4",600,600);

 hframe->Draw();

  // TLine *line2 = new TLine(-7,1,400,1);
  // line->SetLineColor(kGray);
  // line->SetLineWidth(2);
  // line->SetLineStyle(7);
  line->Draw();

 gData_PbPb_b_sys->Draw("e5");
 gData_PbPb_inc_sys->Draw("e5");
 gData_bjt_r->Draw("p");
 gData_inc_r->Draw("p");


  auto lPbpp = new TLegend(0.2,0.2,0.4,0.35);
  lPbpp->AddEntry(gData_bjt_r,"b dijets","P");
  lPbpp->AddEntry(gData_inc_r,"Inclusive dijets","P");
  lPbpp->Draw();




  CMS_lumi(c6, iPeriod, iPos ); 

SavePlot(c6, Pbppratio ? "PbPbppratio" : "PbPbppdiff");

}



  
