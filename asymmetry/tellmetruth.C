#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "TRandom.h"
#include "../corrections/tageffcorrections.h"
#include "../corrections/eclipsecorrections.h"
#include "meanerrors.h"

//for reporting in the end
vector<float> lbinmin, lbinmax;
vector<float> xjdtincmean, xjdtincmeanerror;
vector<float> xjdtbjtmean, xjdtbjtmeanerror;

vector<float> xjmcincmean, xjmcincmeanerror;
vector<float> xjmcincsig2mean, xjmcincsig2meanerror;
vector<float> xjmcbjtSBmean,xjmcbjtSBmeanerror;

vector<float> xjmcbjtmean, xjmcbjtmeanerror;
vector<float> xjdtb12mean, xjdtb12meanerror;
vector<float> xjmcb12mean, xjmcb12meanerror;
vector<float> xjmcb12Signalmean, xjmcb12Signalmeanerror;

TFile *fxjdphi;

TF1 *ftrigEta, *ftrigCent, *ftrigPt;

bool tagLeadingJet = true;
TString dtppjpf = "dtppjpf";
TString dtPbbjt = "dtPbbjt";
TString dtPbjcl = "dtPbjcl";

TString mcppqcd = "mcppqcd";
TString mcppbfa = "mcppbfa";
TString mcPbbfa = "mcPbbfa";
TString mcPbqcd = "mcPbqcd";

bool sampleSubleading = false;
TF1 *fmistagpp = 0, *fmistagfPb1, *fmistagfPb2, *fmistagfPb3;


float NSfracbjt = 1;
bool ppsmearing = false;

void loadmistagsampling()
{
  TFile f("../correctionfiles/BXmistagfunc.root");
  fmistagpp = (TF1 *)f.Get("fpp");
  fmistagfPb1 = (TF1 *)f.Get("fPb1");
  fmistagfPb2 = (TF1 *)f.Get("fPb2");
  fmistagfPb3 = (TF1 *)f.Get("fPb3");
  cout<<"Loaded mistag funcs"<<endl;
}

float getmistagweight(float pt2, int bin)
{
  if (fmistagpp==0) loadmistagsampling();
  if (bin==-1) //pp
    return fmistagpp->Eval(pt2);

  if (bin<20)
    return fmistagfPb1->Eval(pt2);
  if (bin>=20 && bin<60)
    return fmistagfPb2->Eval(pt2);
  if (bin>=60)
    return fmistagfPb3->Eval(pt2);

  return -999;
}

void loadTrigEffCorrections()
{
  TFile *fppFits = new TFile("../correctionfiles/trigEffCorr.root");
   
  ftrigEta = (TF1*) fppFits->Get("fitEta");
  ftrigPt = (TF1*) fppFits->Get("fitPt");
  ftrigCent = (TF1*) fppFits->Get("fitCent");

}

float getTrigcorrection(float jetpt, float jeteta, float bin)
{

  if (!applyTriggerCorr) return 1;

  double centTagEffCorr = 1./ftrigCent->Eval(bin);

  double PtTrigEff = (jetpt > 200) ? ftrigPt->Eval(200):ftrigPt->Eval(jetpt);
  double EtaTrigEff = ftrigEta->Eval(jeteta);

  return centTagEffCorr/PtTrigEff/EtaTrigEff;
}

bool LeadingJetCut(dict &d)
{
  if (tagLeadingJet) return d["discr_csvV1_1"]>csvcut1;
  else return d["discr_csvV1_1"]<0.5;
}

void findtruthpp(float datafraction = 1.)
{
  TFile *fdtpp = new TFile(config.getFileName_djt(dtppjpf));

  TFile *fmcpp = new TFile("/Users/istaslis/Documents/CMS/bjet2015/ntuples/eta1p5/mcppbfaak4PF_djt.root");//config.getFileName_djt(mcppbfa)); //BFA!!!


  seth(10,0,1);
  auto hdtppxJAS = geth("hdtppxJAS","Data b-jets;x_{J}");
  auto hdtINCppxJAS = geth("hdtINCppxJAS","Data Inclusive;x_{J}");
  auto hdt12ppxJAS = geth("hdt12ppxJAS","Data b-jets;x_{J}");

  auto hmcppxJAS = geth("hmcppxJAS","MC b-jets;x_{J}");
  auto hmcppxJASsigBB = geth("hmcppxJASsigBB","MC b-jets;x_{J}");
  auto hmcppxJASsignal2 = geth("hmcppxJASsignal2","MC b-jets;x_{J}");
  auto hmcppqcdxJAS = geth("hmcppqcdxJAS","MC Inclusive;x_{J}");
  auto hmcppqcdxJASsignal2 = geth("hmcppqcdxJASsignal2","MC Inclusive;x_{J}");
  auto hmc12ppxJAS = geth("hmc12ppxJAS","MC b-jets;x_{J}");
  auto hmc12ppxJNS = geth("hmc12ppxJNS","MC b-jetsNS;x_{J}");

  //dphi
  seth(20,0,3.142);
  auto hdphippINCdata = geth("hdphippINCdata","Data Inclusive;#Delta#phi");
  auto hdphippINCmc = geth("hdphippINCmc","MC Inclusive;#Delta#phi");

  auto hdphippBJTdata = geth("hdphippBJTdata","Data b-jets;#Delta#phi");
  auto hdphippBJTmc = geth("hdphippBJTmc","MC b-jets;#Delta#phi");
  auto hdphippBJT12data = geth("hdphippBJT12data","Data b-jets;#Delta#phi");
  auto hdphippBJT12mc = geth("hdphippBJT12mc","MC b-jets;#Delta#phi");

  //rj
  seth(15,100,250);
  auto hrjincdt = geth("hrjincdt","Data;p_{T,1} [GeV];R_{J}");
  auto hrjincmc = geth("hrjincmc","Pythia 6;p_{T,1} [GeV];R_{J}");
  auto hrjbjtdt = geth("hrjbjtdt","b-jets data;p_{T,1} [GeV];R_{J}");
  auto hrjbjtmc = geth("hrjbjtmc","b-jets MC;p_{T,1} [GeV];R_{J}");

  auto hrjincdtden = geth("hrjincdtden");
  auto hrjincmcden = geth("hrjincmcden");
  auto hrjbjtdtden = geth("hrjbjtdtden");
  auto hrjbjtmcden = geth("hrjbjtmcden");

  TString origpt1 = ppsmearing ? "jtptsansjec1" : "jtpt1";
  TString origpt2 = ppsmearing ? "jtptsansjec2" : "jtpt2";

  Fill(fdtpp,[&] (dict &d) {
    if (d["numTagged"]>6) return;
    float w = d["weight"];
    float corr = tageffcorrectionpp(d[origpt1],d["jteta1"],d[origpt2],d["jteta2"]);//(d["jtpt1"],d["jteta1"],d[jtptSL],d[jtetaSL]);

    float wb = w*corr;
    //if the subleading jet is sampled, the weight is increased by sampling weight
    //and subleading jet must be anti-tagged
    if (sampleSubleading) wb*=getmistagweight(d["jtpt2"],-1);
    bool taggedsubleading = sampleSubleading ? d[discr_csvV1_2]<0.5 : d[discr_csvV1_2]>csvcut2;

    if (d["jtpt1"]>pt1cut && LeadingJetCut(d) && d[jtptSL]>pt2cut)
      hdphippBJTdata->Fill(d[dphiSL1],wb);   

    if (d["jtpt1"]>pt1cut && LeadingJetCut(d) && d[jtptSL]>pt2cut && d[dphiSL1]>PI23)
      hdtppxJAS->Fill(d[jtptSL]/d["jtpt1"],wb);

    if (d["jtpt1"]>pt1cut && LeadingJetCut(d) && d["jtpt2"]>pt2cut && taggedsubleading && d["dphi21"]>PI23) {
      hdt12ppxJAS->Fill(d["jtpt2"]/d["jtpt1"],wb);
      hrjbjtdt->Fill(d["jtpt1"],wb);
    }
    if (d["jtpt1"]>pt1cut && LeadingJetCut(d)) hrjbjtdtden->Fill(d["jtpt1"],wb); 

    if (d["jtpt1"]>pt1cut && LeadingJetCut(d) && d["jtpt2"]>pt2cut && taggedsubleading)
      hdphippBJT12data->Fill(d["dphi21"],wb);

    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hdtINCppxJAS->Fill(d["jtpt2"]/d["jtpt1"],w);
      hrjincdt->Fill(d["jtpt1"],w);
    }
     if (d["jtpt1"]>pt1cut) hrjincdtden->Fill(d["jtpt1"],w);

    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut)
      hdphippINCdata->Fill(d["dphi21"],w);


    },datafraction);

  Fill(fmcpp,[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["pthat"]<pthatcut) return;
    if (m["refpt1"]<50) return;
    float w = weight1SLpp(m);
    float wSB = m["weight"]*processweight((int)m["bProdCode"]);
    float corr = tageffcorrectionpp(m["jtpt1"],m["jteta1"],m["jtpt2"],m["jteta2"]);//(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL]);
    float wb = w*corr;
    float w2 = weight12(m);

    //if the subleading jet is sampled, the weight is increased by sampling weight
    //and subleading jet must be anti-tagged
    if (sampleSubleading) wb*=getmistagweight(m["jtpt2"],-1);
    bool taggedsubleading = sampleSubleading ? m[discr_csvV1_2]<0.5 : m[discr_csvV1_2]>csvcut2;

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m[jtptSL]>pt2cut)
      hdphippBJTmc->Fill(m[dphiSL1],wb);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m[jtptSL]>pt2cut && m[dphiSL1]>PI23)
      hmcppxJAS->Fill(m[jtptSL]/m["jtpt1"],wb);

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading)
      hdphippBJT12mc->Fill(m["dphi21"],wb);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]>PI23)
      hmc12ppxJAS->Fill(m["jtpt2"]/m["jtpt1"],wb);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]<PI13)
      hmc12ppxJNS->Fill(m["jtpt2"]/m["jtpt1"],wb);


    if (m["jtpt1"]>pt1cut && m["pairCodeSB1"]==0 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcppxJASsigBB->Fill(m["jtptSB"]/m["jtpt1"],wSB);

    if (m["jtpt1"]>pt1cut && m["pairCodeSignal21"]==0 && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23) {
      hmcppxJASsignal2->Fill(m["jtptSignal2"]/m["jtpt1"],wSB);//weight assumes both are b-jets
      hrjbjtmc->Fill(m["jtpt1"],w2);
    }
    if (m["jtpt1"]>pt1cut && abs(m["refparton_flavorForB1"])==5)
      hrjbjtmcden->Fill(m["jtpt1"],w2);

  },datafraction);

  auto b12sub = (TH1F *)hmc12ppxJAS->Clone("b12sub");
  b12sub->Add(hmc12ppxJNS,-1);
  cout<<"pp AS: "<<hmc12ppxJAS->GetMean()<<"±"<<hmc12ppxJAS->GetMeanError() <<" - "<<b12sub->GetMean()<<"±"<<b12sub->GetMeanError()<<endl;

  auto fmcppqcd = config.getfile_djt(mcppqcd);
  Fill(fmcppqcd,[&] (dict &m) {
    if (m["pthat"]<pthatcut) return;
    if (m["refpt1"]<50) return;
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut && m["dphi21"]>PI23)
      hmcppqcdxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut)
      hdphippINCmc->Fill(m["dphi21"],w);

    if (m["jtpt1"]>pt1cut && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23) {
      hmcppqcdxJASsignal2->Fill(m["jtptSignal2"]/m["jtpt1"],w);
      hrjincmc->Fill(m["jtpt1"],w);
    }
    if (m["jtpt1"]>pt1cut)hrjincmcden->Fill(m["jtpt1"],w);

  },datafraction);

  SetData({hdtppxJAS,hdtINCppxJAS,hdt12ppxJAS,hdphippINCdata,hdphippBJTdata,hdphippBJT12data});
  SetMC({hmcppxJAS,hmcppxJASsigBB,hmcppxJASsignal2,hmcppqcdxJAS,hmc12ppxJAS,hdphippINCmc,hdphippBJTmc,hdphippBJT12mc});
  SetInc({hdtINCppxJAS,hmcppqcdxJAS,hdphippINCmc,hdphippINCdata});
  SetB({hdtppxJAS,hmcppxJAS,hmcppxJASsigBB,hmcppxJASsignal2,hdt12ppxJAS,hmc12ppxJAS,hdphippBJTmc,hdphippBJT12mc,hdphippBJT12data,hdphippBJTdata});

  // hdt12ppxJAS->SetMarkerColor(darkviolet);  hdt12ppxJAS->SetLineColor(darkviolet);
  // hmc12ppxJAS->SetMarkerColor(darkviolet);  hmc12ppxJAS->SetLineColor(darkviolet);
  // hdphippBJT12data->SetMarkerColor(darkviolet);hdphippBJT12data->SetLineColor(darkviolet);
  // hdphippBJT12mc->SetMarkerColor(darkviolet);hdphippBJT12mc->SetLineColor(darkviolet);
  // hmcppxJASsignal2->SetMarkerColor(darkviolet);hmcppxJASsignal2->SetLineColor(darkviolet);

hrjincdt->Divide(hrjincdt,hrjincdtden,1,1);//,"B"
hrjincmc->Divide(hrjincmc,hrjincmcden,1,1);//,"B"
hrjbjtdt->Divide(hrjbjtdt,hrjbjtdtden,1,1);//,"B"
hrjbjtmc->Divide(hrjbjtmc,hrjbjtmcden,1,1);//,"B"

SetInc({hrjincdt,hrjincmc});
SetMC({hrjincmc});SetData({hrjincdt});

aktstring = "";
plotymin=0;plotymax = 1.;
plotoverwritecolors = false;
plotlegendpos = BottomRight;


plotlegendheader = "Inclusive jets pp";
Draw({hrjincdt,hrjincmc});
plotymin=0;plotymax = 0.5;

plotlegendheader = "b-jets pp";
Draw({hrjbjtdt,hrjbjtmc});

plotoverwritecolors = true;


  NormalizeAllHists();
  plotputmean = true;
  plotytitle = "Event fractions";
  plotdivide = false;
  aktstring += Form("R=0.4 |#eta|<%.1f",etacut);
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2b}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";
  plotymax = 0.4;

  plotdiffmax = 0.02;
  // DrawCompare(hdtppxJAS,hmcppxJAS);
  DrawCompare(hdtppxJAS,hmcppxJASsigBB);

  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  DrawCompare(hdtINCppxJAS,hmcppqcdxJAS);
  // DrawCompare(hdtINCppxJAS,hmcppqcdxJASsignal2);
  
  plotthirdline = "";
  plotymax = 0.65;
  plotputmean = false; plotputwidth = true;
  DrawCompare(hdphippBJT12data,hdphippBJT12mc,"#Delta#phi"); 
  plotputmean = true; plotputwidth = false;  
  plotymax = 0.4;
  plotthirdline = "#Delta#phi>2/3#pi";
  //reco
  //DrawCompare(hdt12ppxJAS,hmc12ppxJAS);   
  DrawCompare(hdt12ppxJAS,hmcppxJASsignal2);

  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2b}>%d GeV", (int)pt1cut, (int)pt2cut);

  plotymax = 0.65;
  plotthirdline = "";

  plotputmean = false; plotputwidth = true;
  plotdiffmax = 0.03;
  DrawCompare(hdphippBJTdata,hdphippBJTmc,"#Delta#phi");

  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  DrawCompare(hdphippINCdata,hdphippINCmc,"#Delta#phi");
  plotdiffmax = 9999;
  plotputmean = true; plotputwidth = false;


  lbinmin.push_back(-1);lbinmax.push_back(-1);
  xjdtincmean.push_back(       hdtINCppxJAS->GetMean());
  xjdtincmeanerror.push_back(  hdtINCppxJAS->GetMeanError());
  xjmcincmean.push_back(       hmcppqcdxJAS->GetMean());
  xjmcincmeanerror.push_back(  hmcppqcdxJAS->GetMeanError());
  xjdtbjtmean.push_back(       hdtppxJAS->GetMean());
  xjdtbjtmeanerror.push_back(  hdtppxJAS->GetMeanError());
  xjdtb12mean.push_back(       hdt12ppxJAS->GetMean());
  xjdtb12meanerror.push_back(  hdt12ppxJAS->GetMeanError());
  xjmcbjtmean.push_back(       hmcppxJAS->GetMean());
  xjmcbjtmeanerror.push_back(  hmcppxJAS->GetMeanError());
  xjmcb12mean.push_back(       hmc12ppxJAS->GetMean());
  xjmcb12meanerror.push_back(  hmc12ppxJAS->GetMeanError());
  xjmcb12Signalmean.push_back(       hmc12ppxJAS->GetMean());
  xjmcb12Signalmeanerror.push_back(  hmc12ppxJAS->GetMeanError());
  xjmcbjtSBmean.push_back(         hmcppxJASsigBB->GetMean());//hmcxJASsigBB->GetMean());
  xjmcbjtSBmeanerror.push_back(    hmcppxJASsigBB->GetMeanError());//hmcxJASsigBB->GetMeanError());

  xjmcincsig2mean.push_back(        hmcppqcdxJASsignal2->GetMean());
  xjmcincsig2meanerror.push_back(  hmcppqcdxJASsignal2->GetMeanError());

  
  fxjdphi->cd();


  float w,e;
  auto f1 = fitdphi(hdphippINCdata,w,e);
  auto f2 = fitdphi(hdphippINCmc,w,e);
  auto f3 = fitdphi(hdphippBJT12data,w,e);
  auto f4 = fitdphi(hdphippBJT12mc,w,e);


  hdtINCppxJAS->Write("xJ_data_inc_pp");
  hmcppqcdxJAS->Write("xJ_mc_inc_pp");
  hdt12ppxJAS->Write("xJ_data_bjt_pp");
  hmc12ppxJAS->Write("xJ_mc_bjt_pp");
  hdphippINCdata->Write("dphi_data_inc_pp");
  hdphippINCmc->Write("dphi_mc_inc_pp");
  hdphippBJT12data->Write("dphi_data_bjt_pp");
  hdphippBJT12mc->Write("dphi_mc_bjt_pp");

  f1->Write("fit_dphi_data_inc_pp");
  f2->Write("fit_dphi_mc_inc_pp");
  f3->Write("fit_dphi_data_bjt_pp");
  f4->Write("fit_dphi_mc_bjt_pp");


}


void findtruthPbPb(int binMin, int binMax)
{
  TString dtfname = config.getFileName_djt(dtPbbjt);
  TString mcfname = config.getFileName_djt(mcPbbfa);

  TFile *fdt = new TFile(dtfname);
  TFile *fmc = new TFile(mcfname);
  TFile *fdtinc = new TFile(config.getFileName_djt(dtPbjcl));
  TFile *fmcinc = new TFile(config.getFileName_djt(mcPbqcd));



  buildNamesuffix = TString::Format("_bin_%d_%d",binMin, binMax);
  //  buildTitlesuffix = TString::Format("%d-%d %%",binMin/2, binMax/2);

  //check vertex and centrality first
  seth(40,-15,15);
  auto hdtvz = geth("hdtvz","Data;vz [cm]");
  auto hmcvz = geth("hmcvz","MC;vz [cm]");

  seth(100,0,200);
  auto hdtbin = geth("hdtbin","Data;bin");
  auto hmcbin = geth("hmcbin","MC;bin"); 

  //xJ
  seth(10,0,1);
  auto hdtxJAS = geth("hdtxJAS","Data away-side;x_{J}");
  auto hdtxJNS = geth("hdtxJNS","Data near-side;x_{J}");
  auto hdtxJEars = geth("hdtxJEars","Data ears;x_{J}");
  auto hmcxJAS = geth("hmcxJAS","MC away-side;x_{J}");
  auto hmcxJNS = geth("hmcxJNS","MC near-side;x_{J}");
  auto hmcxJEars = geth("hmcxJEars","MC ears;x_{J}");

  auto hdtxJASSubEars = geth("hdtxJASSubEars","Data b-jets;x_{J}");
  auto hmcxJASSubEars = geth("hmcxJASSubEars","MC b-jets;x_{J}");

  auto hmcxJASsig = geth("hmcxJASsig","SL sig MC;x_{J}");
  auto hmcxJASsigsyslo = geth("hmcxJASsigsyslo","SL sig MC sys lo;x_{J}");
  auto hmcxJASsigsyshi = geth("hmcxJASsigsyshi","SL sig MC sys hi;x_{J}");
  auto hmcxJASsigBB = geth("hmcxJASsigBB","MC b-jets;x_{J}");

  auto hdtxJ12AS = geth("hdtxJ12AS","Data b-jets 12;x_{J}");
  auto hdtxJ12NS = geth("hdtxJ12NS","Data b-jets 12;x_{J}");

  auto hmcxJ12AS = geth("hmcxJ12AS","MC b-jets 12;x_{J}");
  auto hmcxJ12NS = geth("hmcxJ12NS","MC b-jets 12;x_{J}");
  auto hmcxJSignal12AS = geth("hmcxJSignal12AS","MC b-jets 12;x_{J}");


  //inclusive jets
  auto hdtINCxJAS = geth("hdtINCxJAS","INC Data away-side;x_{J}");
  auto hdtINCxJEars = geth("hdtINCxJEars","INC Data ears;x_{J}");
  auto hmcINCxJAS = geth("hmcINCxJAS","INC MC away-side;x_{J}");
  auto hmcINCxJEars = geth("hmcINCxJEars","INC MC ears;x_{J}");

  auto hmcINCxJASsig = geth("hmcINCxJASsig","INC sig MC;x_{J}");

  auto hmcINCxJASsignal2 = geth("hmcINCxJASsignal2","MC Inclusive;x_{J}");


  auto hdtINCxJASSubEars = geth("hdtINCxJASSubEars","Data Inclusive;x_{J}");
  auto hmcINCxJASSubEars = geth("hmcINCxJASSubEars","MC Inclusive;x_{J}");


  //pt2
  seth(10,40,140);
  auto hpt2ASraw = geth("hpt2ASraw","Raw away-side;p_{T,2} [GeV/c]");
  auto hpt2NSraw = geth("hpt2NSraw","Raw near-side;p_{T,2} [GeV/c]");
  auto hpt2AStag = geth("hpt2AStag","+Tageff corrected away-side;p_{T,2} [GeV/c]");
  auto hpt2NStag = geth("hpt2NStag","+Tageff corrected near-side;p_{T,2} [GeV/c]");
  auto hpt2ASecl = geth("hpt2ASecl","+Eclipse corrected away-side;p_{T,2} [GeV/c]");
  auto hpt2NSecl = geth("hpt2NSecl","+Eclipse corrected near-side;p_{T,2} [GeV/c]");
  auto hpt2AStagecl = geth("hpt2AStagecl","+TEC + EC away-side;p_{T,2} [GeV/c]");
  auto hpt2NStagecl = geth("hpt2NStagecl","+TEC + EC near-side;p_{T,2} [GeV/c]");  
  auto hpt2sub = geth("hpt2sub","Subtacted;p_{T,2} [GeV/c]");




  //dphi
  seth(20,0,3.142);
  auto hdphiINCdata = geth("hdphiINCdata","Data Inclusive;#Delta#phi");
  auto hdphiINCall = geth("hdphiINCall","MC Inclusive;#Delta#phi");
  auto hdphiINCsig = geth("hdphiINCsig","MC Inclusive, signal;#Delta#phi");

  auto hdphiBJTdata = geth("hdphiBJTdata","Data b-jets;#Delta#phi");
  auto hdphiBJTall = geth("hdphiBJTall","MC b-jets;#Delta#phi");
  auto hdphiBJTsig = geth("hdphiBJTsig","MC b-jets, signal;#Delta#phi");

  auto hdtdphiBJT12 = geth("hdtdphiBJT12","Data b-jets;#Delta#phi");
  auto hmcdphiBJT12 = geth("hmcdphiBJT12","MC b-jets;#Delta#phi");
  auto hmcdphiBJT12sig = geth("hmcdphiBJT12sig","MC b-jets, signal;#Delta#phi");

  //pair codes
  seth(5,0,5);
  auto hPairCodeQCD = geth("hPairCodeQCD");
  auto hPairCodeBFA = geth("hPairCodeBFA");
  auto hPairCode = geth("hPairCode");

  auto hPairCodeQCD12 = geth("hPairCodeQCD12");
  auto hPairCodeBFA12 = geth("hPairCodeBFA12");
  auto hPairCode12 = geth("hPairCode12");

  //centrality
  seth(10,0,100); //don't forget to divide!!!
  auto hbinconfusion12 = geth("hbinconfusion12","12 analysis;bin;Jet Confusion");
  auto hbinconfusionSL = geth("hbinconfusionSL","SL analysis;bin;Jet Confusion");
  auto hbinSignal = geth("hbinSignal");
  auto hbinSignalFound12 = geth("hbinSignalFound12");
  auto hbinSignalFoundSL = geth("hbinSignalFoundSL");

  auto hbinfakerateSL = geth("hbinfakerateSL","SL analysis;bin;Purity");
  auto hbinSL = geth("hbinSL");
  auto hbinSLisB = geth("hbinSLisB");

  seth(15,100,250);
  auto hrjincdtAS = geth("hrjincdtAS","Data;p_{T,1} [GeV];R_{J}");
  auto hrjincmcAS = geth("hrjincmcAS","Pythia 6 + Hydjet;p_{T,1} [GeV];R_{J}");
  auto hrjbjtdtAS = geth("hrjbjtdtAS","b-jets data;p_{T,1} [GeV];R_{J}");
  auto hrjbjtmcAS = geth("hrjbjtmcAS","b-jets MC;p_{T,1} [GeV];R_{J}");

  auto hrjincdtNS = geth("hrjincdtNS");
  auto hrjincmcNS = geth("hrjincmcNS");
  auto hrjbjtdtNS = geth("hrjbjtdtNS");
  auto hrjbjtmcNS = geth("hrjbjtmcNS");


  auto hrjincdtdenPb = geth("hrjincdtdenPb");
  auto hrjincmcdenPb = geth("hrjincmcdenPb");
  auto hrjbjtdtdenPb = geth("hrjbjtdtdenPb");
  auto hrjbjtmcdenPb = geth("hrjbjtmcdenPb");


  Fill(fdt,[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>=binMax) return;

    float w = m["weight"];
    float corr = tageffcorrectionPbPb(m["jtpt1"],m["jteta1"],m["jtpt2"],m["jteta2"],m["bin"]);
    float trigcorr = getTrigcorrection(m["jtpt1"],m["jteta1"],m["bin"]); //ideally use trigger pt
    float ecorr = eclipseWeightdt(m["jtpt2"],m["bin"]);
    float wb = w*corr*trigcorr;

    if (sampleSubleading) wb*=getmistagweight(m["jtpt2"],m["bin"]);
    bool taggedsubleading = sampleSubleading ? m[discr_csvV1_2]<0.5 : m[discr_csvV1_2]>csvcut2;



    float dphi = m[dphiSL1];
    float deta = abs(m["jteta1"]-m[jtetaSL]);
    hdtvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m[jtptSL]>pt2cut) {

      hdtbin->Fill(m["bin"],wb);
      hdphiBJTdata->Fill(m[dphiSL1],wb);

      if (dphi>PI23)
        hdtxJAS->Fill(m[jtptSL]/m["jtpt1"],wb);
      if (dphi<PI13)
        hdtxJNS->Fill(m[jtptSL]/m["jtpt1"],wb);

      if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
        hdtxJEars->Fill(m[jtptSL]/m["jtpt1"],wb);

    }

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading)
      hdtdphiBJT12->Fill(m["dphi21"],wb); //not eclipse-corrected!

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]>PI23) {
      hdtxJ12AS->Fill(m["jtpt2"]/m["jtpt1"],wb*ecorr);
      hpt2ASraw->Fill(m["jtpt2"],w);
      hpt2AStag->Fill(m["jtpt2"],wb);
      hpt2ASecl->Fill(m["jtpt2"],w*ecorr);
      hpt2AStagecl->Fill(m["jtpt2"],wb*ecorr);
      hrjbjtdtAS->Fill(m["jtpt1"],wb*ecorr);
    }

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]<PI13) {
      hdtxJ12NS->Fill(m["jtpt2"]/m["jtpt1"],wb*ecorr*NSfracbjt);
      hpt2NSraw->Fill(m["jtpt2"],w);
      hpt2NStag->Fill(m["jtpt2"],wb);
      hpt2NSecl->Fill(m["jtpt2"],w*ecorr*NSfracbjt);
      hpt2AStagecl->Fill(m["jtpt2"],wb*ecorr*NSfracbjt);
      hrjbjtdtNS->Fill(m["jtpt1"],wb*ecorr);
    }
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m)) hrjbjtdtdenPb->Fill(m["jtpt1"],wb);

  });



  vector<float> cbins = {0.,20.,60.,200.};
  vector<TString> binnames = {"0-10%", "10-30%", "30-100%"};



  Fill(fmc,[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>=binMax) return;
    if (m["pthat"]<pthatcut) return;
    if (m["event"]==128751 || m["event"]==1551232 || m["event"]==3043866) return;

    //at least one of the two jets must be a b-jet
    if (abs(m["refparton_flavorForB1"])!=5 && abs(m[refparton_flavorForBSL])!=5) return;

    //float w0 = m["weight"];


    float dphi = m[dphiSL1];
    float deta = abs(m["jteta1"]-m[jtetaSL]);

    float w0=weight1SLPbPb(m);
    float corr = tageffcorrectionPbPb(m["jtpt1"],m["jteta1"],m["jtpt2"],m["jteta2"],m["bin"]);//(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL],m["bin"]);
    float w = w0*corr;

    float ecorr = eclipseWeightmc(m["jtpt2"],m["bin"]);

    float w12 = weight12(m);

    float wSB = m["weight"]*processweight((int)m["bProdCode"]);
    float wS2 = m["weight"];
    if (m["pairCodeSignal21"]==0) wS2*=processweight((int)m["bProdCode"]);

    if (sampleSubleading) w*=getmistagweight(m["jtpt2"],m["bin"]);
    bool taggedsubleading = sampleSubleading ? m[discr_csvV1_2]<0.5 : m[discr_csvV1_2]>csvcut2;


    hmcvz->Fill(m["vz"],w);

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading)
      hmcdphiBJT12->Fill(m["dphi21"],w);

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]>PI23) {
      hmcxJ12AS->Fill(m["jtpt2"]/m["jtpt1"],w*ecorr);
      hrjbjtmcAS->Fill(m["jtpt1"],w12*ecorr);
    }
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && taggedsubleading && m["dphi21"]<PI13) {
      hmcxJ12NS->Fill(m["jtpt2"]/m["jtpt1"],w*ecorr*NSfracbjt);
      hrjbjtmcNS->Fill(m["jtpt1"],w12*ecorr);
    }
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m))
      hrjbjtmcdenPb->Fill(m["jtpt1"],w12);
    
    if (m["jtpt1"]>pt1cut && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23 && m["pairCodeSignal21"]==0)
      hmcxJSignal12AS->Fill(m["jtptSignal2"]/m["jtpt1"],wSB);


    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtptSignal2"]>pt2cut && m["discr_csvV1_Signal2"]>csvcut2 && m["dphiSignal21"]>PI23)
      hPairCodeBFA12->Fill(m["pairCodeSignal21"],wS2);

//NOT CORRECTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m) && m["jtptSignal2"]>pt2cut && m[discr_csvV1_Signal2]>csvcut2) { //&& m["dphiSignal21"]>PI23 
      hbinSignal->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==2)
        hbinSignalFound12->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==m[SLord])
        hbinSignalFoundSL->Fill(m["bin"]/2,w);
    }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && LeadingJetCut(m) && m[jtptSL]>pt2cut && m[dphiSL1]>PI23) {
      hbinSL->Fill(m["bin"],w);//wb);
      if (abs(m[refparton_flavorForBSL])==5)
        hbinSLisB->Fill(m["bin"],w);//wb);
    }


    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcxJASsigBB->Fill(m["jtptSB"]/m["jtpt1"],wSB);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m) && m[jtptSL]>pt2cut) {
      hmcbin->Fill(m["bin"],w);//wb);
      hdphiBJTall->Fill(m[dphiSL1],w);//wb);

      if (m[pairCodeSL1]==0)
        hdphiBJTsig->Fill(m[dphiSL1],w);//wb);

      if (m[dphiSL1]>PI23) {
        hmcxJAS->Fill(m[jtptSL]/m["jtpt1"], w);//m[pairCodeSL1]==0 ? wb : wbkg);


        //signal-like pairCode (for MC purity)
        if (IsSignal(m))//m[subidSL]==0)
          hPairCodeBFA->Fill(m[pairCodeSL1],w0);//wb);////////////////////////////////there were w0 here
      }

      if (m[dphiSL1]<PI13)
        hmcxJNS->Fill(m[jtptSL]/m["jtpt1"],w);//wb);


    }

if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m)  && m[jtptSL]>pt2cut 
        && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1)))
        hmcxJEars->Fill(m[jtptSL]/m["jtpt1"],w);

  });

  cout<<"hmcxJ12NS integral/hmcxJ12AS integral = "<<hmcxJ12NS->Integral()/hmcxJ12AS->Integral()<<endl;
  cout<<"hdtxJNS integral/hdtxJAS integral = "<<hdtxJNS->Integral()/hdtxJAS->Integral()<<endl;

  hmcxJ12AS->Add(hmcxJ12NS,-1);


  Fill(fdtinc,[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>=binMax) return;
    if (m["hiHF"]>5500) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);
    float ew = eclipseWeightdt(m["jtpt2"],m["bin"]);
    float nsfraction = NSfracdt(m["bin"]);


    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCdata->Fill(m["dphi21"],w);

      if (dphi>PI23) {
        hdtINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w*ew);
        hrjincdtAS->Fill(m["jtpt1"],w*ew);
      }
      
      // if (dphi<PI13)
      //   hdtINCxJNS->Fill(m["jtpt2"]/m["jtpt1"],w);
    // if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
      if (dphi<PI13) {
        hdtINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w*ew*nsfraction); //USE SIMPLE NEAR SIDE FOR INC JETS
        hrjincdtNS->Fill(m["jtpt1"],w*ew);
      }
    }

    if (m["jtpt1"]>pt1cut)
      hrjincdtdenPb->Fill(m["jtpt1"],w);


  });


  Fill(fmcinc,[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>=binMax) return;
    if (m["pthat"]<pthatcut) return;
    if (m["refpt1"]<50) return;
    if (m["event"]==128751) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    float ew = eclipseWeightmc(m["jtpt2"],m["bin"]);
    float nsfraction = NSfracmc(m["bin"]);


    if (m["jtpt1"]>pt1cut)  hrjincmcdenPb->Fill(m["jtpt1"],w);
    if (m["jtpt1"]>pt1cut && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23)
      hmcINCxJASsignal2->Fill(m["jtptSignal2"]/m["jtpt1"],w);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCall->Fill(m["dphi21"],w);

      if (m["subid2"]==0)
        hdphiINCsig->Fill(m["dphi21"],w);



      if (dphi>PI23) {
          hmcINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w*ew);

          hrjincmcAS->Fill(m["jtpt1"],w*ew);

          if (m["subid2"]==0)
            hmcINCxJASsig->Fill(m["jtpt2"]/m["jtpt1"],w);
      }
      
      // if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
      if (dphi<PI13) {
        hmcINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w*ew*nsfraction); //USE SIMPLE NS for INC jets
        hrjincmcNS->Fill(m["jtpt1"],w*ew);        
      }
    }

    if (m["numTagged"]>6) return;

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m[jtptSL]>pt2cut && m[dphiSL1]>PI23 && IsSignal(m))
        hPairCodeQCD->Fill(m[pairCodeSL1],w);

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["discr_csvV1_Signal2"]>csvcut2 && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23)
        hPairCodeQCD12->Fill(m["pairCodeSL1"],w);

  });


  hPairCode->SetBinContent(1,hPairCodeBFA->GetBinContent(1));
  hPairCode->SetBinContent(2,hPairCodeQCD->GetBinContent(2));
  hPairCode->SetBinContent(3,hPairCodeBFA->GetBinContent(3));
  hPairCode->SetBinContent(4,hPairCodeBFA->GetBinContent(4));
  hPairCode->SetBinContent(5,hPairCodeQCD->GetBinContent(5));

  hPairCode12->SetBinContent(1,hPairCodeBFA12->GetBinContent(1));
  hPairCode12->SetBinContent(2,hPairCodeQCD12->GetBinContent(2));
  hPairCode12->SetBinContent(3,hPairCodeBFA12->GetBinContent(3));
  hPairCode12->SetBinContent(4,hPairCodeBFA12->GetBinContent(4));
  hPairCode12->SetBinContent(5,hPairCodeQCD12->GetBinContent(5));

  //TODO: fix bin handling
  float coef = bkgfractionInNearSide[getbinindex((binMin+binMax)/2)];
  hdtxJASSubEars->Add(hdtxJAS,hdtxJEars,1,-1*coef);
  hmcxJASSubEars->Add(hmcxJAS,hmcxJEars,1,-1*coef);

  cout<<"AS : "<<hmcxJAS->GetMean()<<endl;
  cout<<"NS : "<<hmcxJNS->GetMean()<<endl;
  cout<<"Sub: "<<hmcxJASSubEars->GetMean()<<endl;

  cout<<"   Fraction of data left after subtraction in bin "<<binMin/2<<" - "<<binMax/2<<" = "<<hdtxJASSubEars->Integral()/hdtxJAS->Integral()<<endl;

  hdtINCxJASSubEars->Add(hdtINCxJAS,hdtINCxJEars,1,-1);
  hmcINCxJASSubEars->Add(hmcINCxJAS,hmcINCxJEars,1,-1);


  Draw({hdtxJ12AS,hdtxJ12NS});

  float xjbjt12mean, xjbjt12meanerror;
  getsubmeanerror(hdtxJ12AS,hdtxJ12NS,xjbjt12mean,xjbjt12meanerror);

  hdtxJ12AS->Add(hdtxJ12NS,-1);
  hpt2sub->Add(hpt2AStagecl,hpt2NStagecl,1,-1);


  Draw({hpt2ASraw,hpt2NSraw});
  Draw({hpt2AStag,hpt2NStag});
  Draw({hpt2ASecl,hpt2NSecl});
  Draw({hpt2sub});


  cout<<"COMPARE xJ bjt 12 mean/error!!! bin "<<binMin/2<<" "<<binMax/2<<endl;
  cout<<"hist difference = "<<hdtxJ12AS->GetMean()<<" ± "<<hdtxJ12AS->GetMeanError()<<endl;
  cout<<"mine difference = "<<xjbjt12mean<<" ± "<<xjbjt12meanerror<<endl;

  Draw({hrjbjtmcAS,hrjbjtmcNS,hrjbjtmcdenPb});

  hrjincdtAS->Add(hrjincdtNS,-1);
  hrjincmcAS->Add(hrjincmcNS,-1);
  hrjbjtdtAS->Add(hrjbjtdtNS,-1);
  hrjbjtmcAS->Add(hrjbjtmcNS,-1);

  hrjincdtAS->Divide(hrjincdtAS,hrjincdtdenPb,1,1);//,"B"
  hrjincmcAS->Divide(hrjincmcAS,hrjincmcdenPb,1,1);//,"B"
  hrjbjtdtAS->Divide(hrjbjtdtAS,hrjbjtdtdenPb,1,1);//,"B"
  hrjbjtmcAS->Divide(hrjbjtmcAS,hrjbjtmcdenPb,1,1);//,"B"


  SetInc({hrjincdtAS,hrjincmcAS});
  SetMC({hrjincmcAS}); SetData({hrjincdtAS});

  plotytitle = "";
  plotlegendpos = BottomRight;
  plotymin = 0; plotymax = 1;
  aktstring = "";
  plotsecondline = "";
  plotthirdline = "";
  plotoverwritecolors = false;

  plotputmean = false;
  plotlegendheader = TString::Format("Inclusive jets %d-%d %%",binMin/2, binMax/2);

  Draw({hrjincdtAS,hrjincmcAS});
  
  // Draw({hrjbjtmcAS});
  // Draw({hrjbjtdtAS});

  plotlegendheader = TString::Format("b-jets %d-%d %%",binMin/2, binMax/2);
  Draw({hrjbjtdtAS,hrjbjtmcAS});
  plotlegendheader = "";


  plotoverwritecolors = true; 



  hbinconfusion12->Divide(hbinSignalFound12,hbinSignal,1,1,"B");
  hbinconfusionSL->Divide(hbinSignalFoundSL,hbinSignal,1,1,"B");

  hbinfakerateSL->Divide(hbinSLisB,hbinSL,1,1,"B");

  Normalize({hPairCode});
  Print(hPairCode);

  Normalize({hPairCode12});
  Print(hPairCode12);

  float purity = hPairCode->GetBinContent(1);
  float purity12 = hPairCode12->GetBinContent(1);

  plotylog = false;
  plotdivide = false;
  plotymin = 0.8;
  //plotymax = 0.4;
  aktstring = Form("anti-k_{T} Pu R=0.4 |#eta|<%.1f",etacut);
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2b}>%d GeV", (int)pt1cut, (int)pt2cut);
  
  plotytitle = "Event fractions";

  plottextposx = 10;
  plotlegendpos = TopRight;
  plotputmean = false;
  plotymin = 0.;


// for (int i=0;i<Nb;i++)
  // cout<<"Amount of background in the subtraction "<<binnames[i]<<" : "<<
      // setprecision(2)<<(1-hptNSsig[i]->Integral(0,hptNSsig[i]->GetNbinsX()+1)/(float)hptNS[i]->Integral(0,hptNS[i]->GetNbinsX()+1))*100<<" % "<<endl;;


  plotymin = 0.5;//0.4;
  plotymax = 0.8;//0.8;
  plotlegendpos = BottomRight;
  plottextposx = 0.5;
  plottextposy = 0.79;

  plotputmean = false;


  plotputmean = true;

  plotymin = 0.7;
  plotymax = 1;

  plotthirdline = "#Delta#phi>2/3#pi";

  plottextposx = 0.4;
  plottextposy = 0.65;
   // Draw({hbinconfusion12,hbinconfusionSL});
    plotymin = 0;
     plotlegendpos = None;
   // Draw({hbinfakerateSL});


  plottextposx = 0.55;
  plottextposy = 0.79;


  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity);
  plotlegendpos = TopLeft;

  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcINCxJASSubEars});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsig});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsig});

  plotputmean = true;
  plotputwidth = false;
  plotymax = 9999;

  // Draw({hmcINCxJAS,hmcINCxJEars,hmcINCxJASSubEars,hmcINCxJASsig});
  // Draw({hmcxJAS,hmcxJEars,hmcxJASSubEars,hmcxJASsig});

  // DrawCompare(hmcINCxJAS,hmcINCxJEars);
  // DrawCompare(hmcINCxJASSubEars,hmcINCxJASsig);
  // DrawCompare(hmcxJAS,hmcxJEars);
  // DrawCompare(hmcxJASSubEars,hmcxJASsig);

  
  SetMC({hmcdphiBJT12,hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdtdphiBJT12,hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hmcdphiBJT12,hdtdphiBJT12,hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtxJ12AS,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcxJASsigBB,hmcxJSignal12AS,hmcINCxJASSubEars,hmcINCxJASsignal2});

  SetB({hdtxJASSubEars,hdtxJ12AS,hmcxJASSubEars,hmcxJASsigBB,hmcxJSignal12AS,hmcxJASsig,hmcxJASsigsyshi,hmcxJASsigsyslo});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsignal2,hmcINCxJASsig});


  //pt2 distribution shouldn't be normalized
  Normalize({hdphiBJTdata,hdphiBJTall,hdtdphiBJT12,hmcdphiBJT12,hdtxJASSubEars, hmcxJASsigBB,
              hdtxJ12AS,hmcxJ12AS,hdtxJ12AS,hmcxJSignal12AS,hdtINCxJASSubEars,hmcINCxJASsignal2,hdphiINCdata,hdphiINCall,hPairCode});

  plotputmean = true;
  plotymax = 0.4;
    // Draw({hmcxJASsig,hmcxJASsigsyshi,hmcxJASsigsyslo});

  plotputmean = false;
  plotputwidth = true;


  plotthirdline = TString::Format("%d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity);
  plotdiffmax = 0.055;
  plotymax = 0.5;
  DrawCompare(hdphiBJTdata,hdphiBJTall,"#Delta#phi");

  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = TString::Format("%d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity12);
  DrawCompare(hdtdphiBJT12,hmcdphiBJT12,"#Delta#phi");
  
  plotdiffmax = 9999;  
  // DrawCompare(hdphiBJTdata,hdphiBJTsig);

  // Draw({hdphiINCdata,hdphiINCall,hdphiINCsig});
  // Draw({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});
  plotputmean = true;
  plotputwidth = false;

  // plotputmean = true;
  // plotymax = 0.26;

  // Draw({hmcINCxJASsig,hmcINCxJEars,hmcINCxJAS,hmcINCxJASSubEars});
  // Draw({hmcxJASsig,hmcxJEars,hmcxJAS,hmcxJASSubEars});

  plotymax = 0.4;
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2b}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity);
  plotdiffmax = 0.15;
  //important! but not used in the result plots
  // DrawCompare(hdtxJASSubEars, hmcxJASSubEars);
  DrawCompare(hdtxJASSubEars, hmcxJASsigBB);

  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity12);

  //reco
  DrawCompare(hdtxJ12AS,hmcxJ12AS);
  //signal
  DrawCompare(hdtxJ12AS,hmcxJSignal12AS);



  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %%",binMin/2, binMax/2);


  // DrawCompare(hdtINCxJASSubEars,hmcINCxJASSubEars);
  DrawCompare(hdtINCxJASSubEars,hmcINCxJASsignal2);
plotdiffmax = 9999;

  plotputmean = false;
  plotputwidth = true;
  plotthirdline = TString::Format("%d-%d %%",binMin/2, binMax/2);
  plotdiffmax = 0.07;
      plotymax = 0.65;
  DrawCompare(hdphiINCdata,hdphiINCall,"#Delta#phi");

  plotputmean = true;
  plotputwidth = false;

  // DrawCompare(hdphiINCdata,hdphiINCsig);

  plotputmean = false;
  plotymax = 1;
  Draw({hPairCode});

  plotytitle = "";

  // Draw({hdtxJAS,hmcxJAS});
  // Draw({hdtxJNS,hmcxJNS});
  // Draw({hdtxJASSub,hmcxJASSub});
  // Draw({hdtxJASSubEars,hmcxJASSubEars,hdtINCxJASSub});
  // Draw({hdtxJASSubEars,hmcxJASSubEars});
  // Draw({hdtxJASSubEars});
  // Draw({hdtINCxJASSub});
  // Draw({hdtxJASSubEars,hdtINCxJASSub});


  lbinmin.push_back(binMin);lbinmax.push_back(binMax);
  xjdtincmean.push_back(       hdtINCxJASSubEars->GetMean());
  xjdtincmeanerror.push_back(  hdtINCxJASSubEars->GetMeanError());
  xjmcincmean.push_back(       hmcINCxJASSubEars->GetMean());
  xjmcincmeanerror.push_back(  hmcINCxJASSubEars->GetMeanError());
  xjdtbjtmean.push_back(       hdtxJASSubEars->GetMean());
  xjdtbjtmeanerror.push_back(  hdtxJASSubEars->GetMeanError());
  xjdtb12mean.push_back(       hdtxJ12AS->GetMean());
  xjdtb12meanerror.push_back(  hdtxJ12AS->GetMeanError());
  xjmcb12mean.push_back(       hmcxJ12AS->GetMean());
  xjmcb12meanerror.push_back(  hmcxJ12AS->GetMeanError());
  xjmcb12Signalmean.push_back(       hmcxJSignal12AS->GetMean());
  xjmcb12Signalmeanerror.push_back(  hmcxJSignal12AS->GetMeanError());
  xjmcbjtmean.push_back(       hmcxJASSubEars->GetMean());
  xjmcbjtmeanerror.push_back(  hmcxJASSubEars->GetMeanError());
  xjmcbjtSBmean.push_back(         hmcxJASsigBB->GetMean());
  xjmcbjtSBmeanerror.push_back(    hmcxJASsigBB->GetMeanError());
  xjmcincsig2mean.push_back(       hmcINCxJASsignal2->GetMean());
  xjmcincsig2meanerror.push_back(  hmcINCxJASsignal2->GetMeanError());
  

  fxjdphi->cd();

  float w,e;
  auto f1 = fitdphi(hdphiINCdata,w,e);
  auto f2 = fitdphi(hdphiINCall,w,e);
  auto f3 = fitdphi(hdtdphiBJT12,w,e);
  auto f4 = fitdphi(hmcdphiBJT12,w,e);

  hdtINCxJASSubEars->Write(Form("xJ_data_inc_%d_%d",binMin/2,binMax/2));
  hmcINCxJASSubEars->Write(Form("xJ_mc_inc_%d_%d",binMin/2,binMax/2));
  hdphiINCdata->Write(Form("dphi_data_inc_%d_%d",binMin/2,binMax/2));
  f1->Write(Form("fit_dphi_data_inc_%d_%d",binMin/2,binMax/2));
  hdphiINCall->Write(Form("dphi_mc_inc_%d_%d",binMin/2,binMax/2));
  f2->Write(Form("fit_dphi_mc_inc_%d_%d",binMin/2,binMax/2));


  hdtxJ12AS->Write(Form("xJ_data_bjt_%d_%d",binMin/2,binMax/2));
  hmcxJSignal12AS->Write(Form("xJ_mc_bjt_%d_%d",binMin/2,binMax/2));
  hdtdphiBJT12->Write(Form("dphi_data_bjt_%d_%d",binMin/2,binMax/2));
  f3->Write(Form("fit_dphi_data_bjt_%d_%d",binMin/2,binMax/2));
  hmcdphiBJT12->Write(Form("dphi_mc_bjt_%d_%d",binMin/2,binMax/2));
  f4->Write(Form("fit_dphi_mc_bjt_%d_%d",binMin/2,binMax/2));

}

void tellmetruth(TString name = "", int eclmode = 0, int bkgsubtractionmode = 0, bool applytagg = true, bool tagLJ = true, bool tagSJ = true, float csvmode1 = 0.9, float csvmode2 = 0.9, int smearingmode = 0)
{
  name = "results_"+name;
  macro m(name);

  fxjdphi = new TFile(plotfoldername+"/xJdphi.root","recreate");


  bool applytriggercorr = true;
  bool tagSL = true;

  loadTrigEffCorrections();
  inittageffcorr(csvmode1, csvmode2);

  eclipsemode = eclmode;
  bkgsubmode = bkgsubtractionmode;
  NSfracbjt =  1.;//NSfractionbjt;

  tagLeadingJet = tagLJ;

  //tagSL == !mockSL
  dtppjpf = tagSL ? "dtppjpf" : "dXppjpf";
  dtPbbjt = tagSL ? "dtPbbjt" : "dXPbbjt";
  dtPbjcl = tagSL ? "dtPbjcl" : "dXPbj60";

 //mcppqcd
 // mcPbbfa
  if (!tagLJ && !tagSL) {    
    mcppqcd = "mXppqcd";
    mcPbqcd = "mXPbqcd";
    mcPbbfa = "mXPbqcd";
  }


  //in case LJ is not tagged - use inc jet ntuple
  if (!tagLJ && tagSL)  dtPbbjt = "dtPbjcl";
  if (!tagLJ && !tagSL) dtPbbjt = "dXPbj60";

  applyTriggerCorr = applytriggercorr;
  applyCorrection = applytagg;
  sampleSubleading = !tagSJ;

  //USE INCLUSIVE JET NTUPLE!
  // dtPbbjt = "dtPbjcl";
  // applyTriggerCorr=false;

  vector<TString> dtppfiles = {"dtppjpf","dtp1jpf","dtp2jpf","dtp3jpf"};
  vector<TString> mcppfiles = {"mcppqcd","mcp1qcd","mcp2qcd","mcp3qcd"};
  vector<TString> mcppbfafiles = {"mcppbfa","mcp1bfa","mcp2bfa","mcp3bfa"};
  ppsmearing = smearingmode!=0;
  dtppjpf = dtppfiles[smearingmode];
  mcppqcd = mcppfiles[smearingmode];
  mcppbfa = mcppbfafiles[smearingmode];
  findtruthpp(); //fraction of data to process

  findtruthPbPb(60,200);
  findtruthPbPb(20,60); 
  findtruthPbPb(0,20);
  
  //findtruthPbPb(0,200);
  

  cout<<"Bin \t\tInc.Data   \tInc.MC    \tSL.Data   \tSL.MC   \t12.Data   \t12.MC"<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    cout<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t"<<xjdtincmean[i]<<"\t"<<xjdtincmeanerror[i]<<
                                          "\t"<<xjmcincmean[i]<<"\t"<<xjmcincmeanerror[i]<<
                                          "\t"<<xjdtbjtmean[i]<<"\t"<<xjdtbjtmeanerror[i]<<
                                          "\t"<<xjmcbjtmean[i]<<"\t"<<xjmcbjtmeanerror[i]<<
                                          "\t"<<xjdtb12mean[i]<<"\t"<<xjdtb12meanerror[i]<<
                                          "\t"<<xjmcb12mean[i]<<"\t"<<xjmcb12meanerror[i]<<endl;


  std::ofstream ofs (plotfoldername+"/results.txt", std::ofstream::out);

  ofs<<"Bin \t\tInc.Data   \tInc.MC    \tSL.Data   \tSL.MC   \t12.Data   \t12.MC"<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    ofs<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t"<<xjdtincmean[i]<<"\t"<<xjdtincmeanerror[i]<<
                                          "\t\t"<<xjmcincmean[i]<<"\t"<<xjmcincmeanerror[i]<<
                                          "\t\t"<<xjdtbjtmean[i]<<"\t"<<xjdtbjtmeanerror[i]<<
                                          "\t\t"<<xjmcbjtmean[i]<<"\t"<<xjmcbjtmeanerror[i]<<
                                          "\t\t"<<xjdtb12mean[i]<<"\t"<<xjdtb12meanerror[i]<<
                                          "\t\t"<<xjmcb12mean[i]<<"\t"<<xjmcb12meanerror[i]<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    ofs<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t";
  ofs<<endl;
  for (unsigned i=0;i<xjdtbjtmean.size();i++)
    ofs<<setprecision(3)<<xjdtbjtmean[i]<<"\t";
  ofs<<endl;
    
  ofs.close();


  std::ofstream csv (plotfoldername+"/results.csv", std::ofstream::out);
  csv<<"Bin,Inc.Data, Error ,Inc.MC, Error, 12.Data, Error, 12.MC, Error"<<endl;
  for (unsigned i=0;i<lbinmin.size();i++) {
    if (lbinmin[i]==-1) csv<<"pp";
    else csv<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2;
    csv<<setprecision(3)<<","<<xjdtincmean[i]<<","<<xjdtincmeanerror[i]
                        <<","<<xjmcincmean[i]<<","<<xjmcincmeanerror[i]
                     // <<","<<xjdtbjtmean[i]<<","<<xjdtbjtmeanerror[i]
                     // <<","<<xjmcbjtmean[i]<<","<<xjmcbjtmeanerror[i]
                        <<","<<xjdtb12mean[i]<<","<<xjdtb12meanerror[i]
                        <<","<<xjmcb12mean[i]<<","<<xjmcb12meanerror[i]
                        <<endl;
  }
   
  csv.close();



  map<TString, float> res;
  for (unsigned i=0;i<lbinmin.size();i++) {
    const char* end = lbinmin[i]==-1 ? "pp" : Form("%d%d",(int)lbinmin[i],(int)lbinmax[i]);
    res[Form("xj_data_inc_mean%s",end)]      = xjdtincmean[i];
    res[Form("xj_data_inc_meanerror%s",end)] = xjdtincmeanerror[i];

    res[Form("xj_mc_inc_mean%s",end)]      = xjmcincmean[i];
    res[Form("xj_mc_inc_meanerror%s",end)] = xjmcincmeanerror[i];

    res[Form("xj_mc_sig2_inc_mean%s",end)]      = xjmcincsig2mean[i];
    res[Form("xj_mc_sig2_inc_meanerror%s",end)] = xjmcincsig2meanerror[i];

    res[Form("xj_data_bjt_mean%s",end)]      = xjdtbjtmean[i];
    res[Form("xj_data_bjt_meanerror%s",end)] = xjdtbjtmeanerror[i];

    res[Form("xj_data_b12_mean%s",end)]      = xjdtb12mean[i];
    res[Form("xj_data_b12_meanerror%s",end)] = xjdtb12meanerror[i];

    res[Form("xj_mc_b12_mean%s",end)]      = xjmcb12mean[i];
    res[Form("xj_mc_b12_meanerror%s",end)] = xjmcb12meanerror[i];

    res[Form("xj_mc_b12Signal_mean%s",end)]      = xjmcb12Signalmean[i];
    res[Form("xj_mc_b12Signal_meanerror%s",end)] = xjmcb12Signalmeanerror[i];

    res[Form("xj_mc_bjt_mean%s",end)]      = xjmcbjtmean[i];
    res[Form("xj_mc_bjt_meanerror%s",end)] = xjmcbjtmeanerror[i];

    res[Form("xj_mc_bjtSB_mean%s",end)]      = xjmcbjtSBmean[i];
    res[Form("xj_mc_bjtSB_meanerror%s",end)] = xjmcbjtSBmeanerror[i];

  }

  WriteToFile(plotfoldername+"/results.root",res);

  fxjdphi->Close();
  //moneyplot("moneyplot"+suffix);
}
