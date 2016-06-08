#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "TRandom.h"
#include "../corrections/tageffcorrections.h"
#include "../corrections/eclipseclosure.C"
#include "moneyplot.C"

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



TF1 *ftrigEta, *ftrigCent, *ftrigPt;

bool tagLeadingJet = true;
TString dtppjpf = "dtppjpf";
TString dtPbbjt = "dtPbbjt";
TString dtPbj60 = "dtPbj60";

TString mcppqcd = "mcppqcd";
TString mcPbbfa = "mcPbbfa";
TString mcPbqcd = "mcPbqcd";

float etacut = 1.5;

void loadTrigEffCorrections()
{
  TFile *fppFits = new TFile("/Users/istaslis/Documents/CMS/bjet2015/ntuples/trigEffCorr.root");
   
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
  if (tagLeadingJet) return d["discr_csvV1_1"]>0.9;
  else return d["discr_csvV1_1"]<0.5;
}

void findtruthpp(float datafraction = 1.)
{
  TFile *fdtpp = new TFile(config.getFileName_djt(dtppjpf));

  seth(10,0,1);
  auto hdtppxJAS = geth("hdtppxJAS","Data b-jets;x_{J}");
  auto hdtINCppxJAS = geth("hdtINCppxJAS","Data Inclusive;x_{J}");
  auto hdt12ppxJAS = geth("hdt12ppxJAS","Data b-jets;x_{J}");

  auto hmcppxJAS = geth("hmcppxJAS","MC b-jets;x_{J}");
  auto hmcppxJASsigBB = geth("hmcppxJASsigBB","MC b-jets;x_{J}");
  auto hmcppqcdxJAS = geth("hmcppqcdxJAS","MC Inclusive;x_{J}");
  auto hmcppqcdxJASsignal2 = geth("hmcppqcdxJASsignal2","MC Inclusive;x_{J}");
  auto hmc12ppxJAS = geth("hmc12ppxJAS","MC b-jets;x_{J}");

  //dphi
  seth(20,0,3.142);
  auto hdphippINCdata = geth("hdphippINCdata","Data Inclusive;#Delta#phi");
  auto hdphippINCmc = geth("hdphippINCmc","MC Inclusive;#Delta#phi");

  auto hdphippBJTdata = geth("hdphippBJTdata","Data b-jets;#Delta#phi");
  auto hdphippBJTmc = geth("hdphippBJTmc","MC b-jets;#Delta#phi");
  auto hdphippBJT12data = geth("hdphippBJT12data","Data b-jets;#Delta#phi");
  auto hdphippBJT12mc = geth("hdphippBJT12mc","MC b-jets;#Delta#phi");


  Fill(fdtpp,{"weight","jtpt1",discr_csvV1_1,jtptSL,dphiSL1,"jtpt2","dphi21",discr_csvV1_2,"jteta1",jtetaSL,"numTagged"},[&] (dict &_) {
    if (_["numTagged"]>6) return;
    float w = _["weight"];
    float corr = getppcorrection(_["jtpt1"],_["jteta1"],_[jtptSL],_[jtetaSL]);
    float wb = w*corr;

    if (_["jtpt1"]>pt1cut && LeadingJetCut(_) && _[jtptSL]>pt2cut)
      hdphippBJTdata->Fill(_[dphiSL1],wb);   

    if (_["jtpt1"]>pt1cut && LeadingJetCut(_) && _[jtptSL]>pt2cut && _[dphiSL1]>PI23)
      hdtppxJAS->Fill(_[jtptSL]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && LeadingJetCut(_) && _["jtpt2"]>pt2cut && _[discr_csvV1_2]>0.9 && _["dphi21"]>PI23)
      hdt12ppxJAS->Fill(_["jtpt2"]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && LeadingJetCut(_) && _["jtpt2"]>pt2cut && _[discr_csvV1_2]>0.9)
      hdphippBJT12data->Fill(_["dphi21"],w);

    if (_["jtpt1"]>pt1cut && _["jtpt2"]>pt2cut && _["dphi21"]>PI23)
      hdtINCppxJAS->Fill(_["jtpt2"]/_["jtpt1"],w);
    if (_["jtpt1"]>pt1cut && _["jtpt2"]>pt2cut)
      hdphippINCdata->Fill(_["dphi21"],w);


    },datafraction);

  TFile *fmcpp = new TFile(config.getFileName_djt("mcppbfa"));//("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");
  Fill(fmcpp,{"pthat","numTagged","weight","jtpt1","refpt1","jteta1",discr_csvV1_1,jtptSL,dphiSL1,"bProdCode","jtpt2", discr_csvV1_2, "dphi21",jtetaSL,pairCodeSL1,
              "refparton_flavorForB1","jtptSB","dphiSB1"},[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["pthat"]<pthatcut) return;
    float w = weight1SLpp(m);//m["weight"]*processWeights[(int)m["bProdCode"]];
    float wSB = m["weight"]*processweight((int)m["bProdCode"]);
    float corr = getppcorrection(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL]);
    float wb = w*corr;

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m) && m[jtptSL]>pt2cut)
      hdphippBJTmc->Fill(m[dphiSL1],wb);
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m) && m[jtptSL]>pt2cut && m[dphiSL1]>PI23)
      hmcppxJAS->Fill(m[jtptSL]/m["jtpt1"],wb);

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m[discr_csvV1_2]>0.9)
      hdphippBJT12mc->Fill(m["dphi21"],wb);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m[discr_csvV1_2]>0.9 && m["dphi21"]>PI23)
      hmc12ppxJAS->Fill(m["jtpt2"]/m["jtpt1"],wb);


    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcppxJASsigBB->Fill(m["jtptSB"]/m["jtpt1"],wSB);

  });

  TFile *fmcppqcd = new TFile(config.getFileName_djt(mcppqcd));//("/data_CMS/cms/lisniak/bjet2015/mcppqcdak4PF_djt.root");


  Fill(fmcppqcd,{"pthat","weight","jtpt1","refpt1","jtpt2","dphi21","jtptSignal2","dphiSignal21"},[&] (dict &m) {
    if (m["pthat"]<pthatcut) return;
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut && m["dphi21"]>PI23)
      hmcppqcdxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut)
      hdphippINCmc->Fill(m["dphi21"],w);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23)
      hmcppqcdxJASsignal2->Fill(m["jtptSignal2"]/m["jtpt1"],w);

  });

  SetData({hdtppxJAS,hdtINCppxJAS,hdt12ppxJAS,hdphippINCdata,hdphippBJTdata,hdphippBJT12data});
  SetMC({hmcppxJAS,hmcppxJASsigBB,hmcppqcdxJAS,hmc12ppxJAS,hdphippINCmc,hdphippBJTmc,hdphippBJT12mc});
  SetInc({hdtINCppxJAS,hmcppqcdxJAS,hdphippINCmc,hdphippINCdata});
  SetB({hdtppxJAS,hmcppxJAS,hmcppxJASsigBB,hdt12ppxJAS,hmc12ppxJAS,hdphippBJTmc,hdphippBJT12mc,hdphippBJTdata});

  hdt12ppxJAS->SetMarkerColor(darkviolet);  hdt12ppxJAS->SetLineColor(darkviolet);
  hmc12ppxJAS->SetMarkerColor(darkviolet);  hmc12ppxJAS->SetLineColor(darkviolet);
  hdphippBJT12data->SetMarkerColor(darkviolet);hdphippBJT12data->SetLineColor(darkviolet);
  hdphippBJT12mc->SetMarkerColor(darkviolet);hdphippBJT12mc->SetLineColor(darkviolet);

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
  DrawCompare(hdt12ppxJAS,hmc12ppxJAS);   


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

  
  // TFile *f = new TFile("xJdphi_pp.root","recreate");
  // hdtINCppxJAS->Write("xJ_data_inc_pp");
  // hmcppqcdxJAS->Write("xJ_mc_inc_pp");
  // f->Close();

}

void findtruthPbPb(int binMin, int binMax)
{
  TString dtfname = config.getFileName_djt(dtPbbjt);
  TString mcfname = config.getFileName_djt(mcPbbfa);

  TFile *fdt = new TFile(dtfname);
  TFile *fmc = new TFile(mcfname);
  TFile *fdtinc = new TFile(config.getFileName_djt(dtPbj60));
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

  auto hdtxJ12AS = geth("hdtxJ12AS","Data 12 away-side;x_{J}");
  auto hdtxJ12NS = geth("hdtxJ12NS","Data 12 near-side;x_{J}");
  auto hmcxJ12AS = geth("hmcxJ12AS","MC 12 away-side;x_{J}");
  auto hmcxJ12NS = geth("hmcxJ12NS","MC 12 near-side;x_{J}");
  auto hmcxJSignal12AS = geth("hmcxJSignal12AS","MC 12 Signal;x_{J}");


  //inclusive jets
  auto hdtINCxJAS = geth("hdtINCxJAS","INC Data away-side;x_{J}");
  auto hdtINCxJEars = geth("hdtINCxJEars","INC Data ears;x_{J}");
  auto hmcINCxJAS = geth("hmcINCxJAS","INC MC away-side;x_{J}");
  auto hmcINCxJEars = geth("hmcINCxJEars","INC MC ears;x_{J}");

  auto hmcINCxJASsig = geth("hmcINCxJASsig","INC sig MC;x_{J}");

  auto hmcINCxJASsignal2 = geth("hmcINCxJASsignal2","MC Inclusive;x_{J}");


  auto hdtINCxJASSubEars = geth("hdtINCxJASSubEars","Data Inclusive;x_{J}");
  auto hmcINCxJASSubEars = geth("hmcINCxJASSubEars","MC Inclusive;x_{J}");




  //dphi
  seth(20,0,3.142);
  auto hdphiINCdata = geth("hdphiINCdata","Data Inclusive;#Delta#phi");
  auto hdphiINCall = geth("hdphiINCall","MC Inclusive;#Delta#phi");
  auto hdphiINCsig = geth("hdphiINCsig","MC Inclusive, signal;#Delta#phi");
  auto hdphiBJTdata = geth("hdphiBJTdata","Data b-jets;#Delta#phi");
  auto hdphiBJTall = geth("hdphiBJTall","MC b-jets;#Delta#phi");
  auto hdphiBJTsig = geth("hdphiBJTsig","MC b-jets, signal;#Delta#phi");



  //pair codes
  seth(5,0,5);
  auto hPairCodeQCD = geth("hPairCodeQCD");
  auto hPairCodeBFA = geth("hPairCodeBFA");
  auto hPairCode = geth("hPairCode");

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

  Fill(fdt,{"weight","vz","bin","jtpt1",discr_csvV1_1,jtptSL,dphiSL1,"jteta1",jtetaSL,"numTagged","jtpt2","dphi21","discr_csvV1_2"},[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL],m["bin"]);
    float trigcorr = getTrigcorrection(m["jtpt1"],m["jteta1"],m["bin"]); //WRONG! use trigger pt
    float ecorr = eclipseWeightdt(m["jtpt2"],m["bin"]);
    float wb = w*corr*trigcorr;
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

    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]>PI23)
      hdtxJ12AS->Fill(m["jtpt2"]/m["jtpt1"],wb*ecorr);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]<PI13)
      hdtxJ12NS->Fill(m["jtpt2"]/m["jtpt1"],wb*ecorr);
    

  });

  vector<float> cbins = {0.,20.,60.,200.};
  vector<TString> binnames = {"0-10%", "10-30%", "30-100%"};



  Fill(fmc,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1",discr_csvV1_1,jtptSL,dphiSL1,"jteta1",jtetaSL,subidSL,refptSL,pairCodeSL1,
            "jtptSignal2",discr_csvV1_Signal2,"Signal2ord","SLord","dphiSignal21",refparton_flavorForBSL,"refparton_flavorForB1","numTagged","jtptSB","dphiSB1",
            "jtpt2","pairCodeSignal21","discr_csvV1_2","dphi21"},[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["pthat"]<pthatcut) return;

    //at least one of the two jets must be a b-jet
    if (abs(m["refparton_flavorForB1"])!=5 && abs(m[refparton_flavorForBSL])!=5) return;

    //float w0 = m["weight"];


    float dphi = m[dphiSL1];
    float deta = abs(m["jteta1"]-m[jtetaSL]);

    float w0=weight1SLPbPb(m);//w0*processWeights[(int)m["bProdCode"]];

    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL],m["bin"]);
    float w = w0*corr;

    float ecorr = eclipseWeightmc(m["jtpt2"],m["bin"]);

    float wSB = m["weight"]*processweight((int)m["bProdCode"]);

    //float wbkg = w0*corr;//*(m["pthat"]>80);

    hmcvz->Fill(m["vz"],w);



    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]>PI23)
      hmcxJ12AS->Fill(m["jtpt2"]/m["jtpt1"],w*ecorr);
    if (m["jtpt1"]>pt1cut && LeadingJetCut(m) && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]<PI13)
      hmcxJ12NS->Fill(m["jtpt2"]/m["jtpt1"],w*ecorr);
    
    if (m["jtpt1"]>pt1cut && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23 && m["pairCodeSignal21"]==0)
      hmcxJSignal12AS->Fill(m["jtptSignal2"]/m["jtpt1"],wSB);


//NOT CORRECTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m) && m["jtptSignal2"]>pt2cut && m[discr_csvV1_Signal2]>0.9) { //&& m["dphiSignal21"]>PI23 
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

// // here I subtract background with TRUE leading jet!!!

// float LJeff = 0.35;

if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && LeadingJetCut(m)  && m[jtptSL]>pt2cut //&& abs(m["refparton_flavorForB1"])==5
        && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1)))
        hmcxJEars->Fill(m[jtptSL]/m["jtpt1"],w);//LFeff*w//LJeff*(m[pairCodeSL1]==0 ? wb : wbkg));

  });

  hmcxJ12AS->Add(hmcxJ12NS,-1);


  Fill(fdtinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2","hiHF"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["hiHF"]>5500) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);
    float ew = eclipseWeightdt(m["jtpt2"],m["bin"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCdata->Fill(m["dphi21"],w);

      if (dphi>PI23)
        hdtINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w*ew);
      // if (dphi<PI13)
      //   hdtINCxJNS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
        hdtINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w*ew);
    }


  });

  Fill(fmcinc,{"weight","jtpt1","refpt1","jtpt2","dphi21","bin","jteta1","jteta2",pairCodeSL1,discr_csvV1_1,jtptSL,dphiSL1,"pthat",subidSL,refptSL,"subid2","numTagged","jtptSignal2","dphiSignal21"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["pthat"]<pthatcut) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    float ew = eclipseWeightmc(m["jtpt2"],m["bin"]);




    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["jtptSignal2"]>pt2cut && m["dphiSignal21"]>PI23)
      hmcINCxJASsignal2->Fill(m["jtptSignal2"]/m["jtpt1"],w);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["jtpt2"]>pt2cut) {

      hdphiINCall->Fill(m["dphi21"],w);

      if (m["subid2"]==0)
        hdphiINCsig->Fill(m["dphi21"],w);

      if (dphi>PI23) {
          hmcINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w*ew);

          if (m["subid2"]==0)
            hmcINCxJASsig->Fill(m["jtpt2"]/m["jtpt1"],w);
      }
      
      if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
          hmcINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w*ew);
    }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && LeadingJetCut(m) && m[jtptSL]>pt2cut && m[dphiSL1]>PI23) {
      if (IsSignal(m) && m["numTagged"]<=6)//m[pairCodeSL1]<4 && 
        hPairCodeQCD->Fill(m[pairCodeSL1],w);
    }

  });


  hPairCode->SetBinContent(1,hPairCodeBFA->GetBinContent(1));
  hPairCode->SetBinContent(2,hPairCodeQCD->GetBinContent(2));
  hPairCode->SetBinContent(3,hPairCodeBFA->GetBinContent(3));
  hPairCode->SetBinContent(4,hPairCodeBFA->GetBinContent(4));
  hPairCode->SetBinContent(5,hPairCodeQCD->GetBinContent(5));


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


  // Draw({hdtxJ12NS,hdtxJ12AS});

  plotputmean = true;
  hdtxJ12AS->Add(hdtxJ12NS,-1);




  hbinconfusion12->Divide(hbinSignalFound12,hbinSignal,1,1,"B");
  hbinconfusionSL->Divide(hbinSignalFoundSL,hbinSignal,1,1,"B");

  hbinfakerateSL->Divide(hbinSLisB,hbinSL,1,1,"B");

  Normalize({hPairCode});
  Print(hPairCode);
  float purity = hPairCode->GetBinContent(1);

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
  Draw({hbinconfusion12,hbinconfusionSL});
    plotymin = 0;
     plotlegendpos = None;
  Draw({hbinfakerateSL});


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




  
  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcxJASsigBB,hmcINCxJASSubEars,hmcINCxJASsignal2});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsigBB,hmcxJASsig,hmcxJASsigsyshi,hmcxJASsigsyslo});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsignal2,hmcINCxJASsig});


  Normalize({hdtxJ12AS});
  Draw({hdtxJ12AS});


  NormalizeAllHists();

  plotputmean = true;
  plotymax = 0.4;
    // Draw({hmcxJASsig,hmcxJASsigsyshi,hmcxJASsigsyslo});

  plotputmean = false;
  plotputwidth = true;


  plotthirdline = TString::Format("%d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity);
  plotdiffmax = 0.055;
    plotymax = 0.5;
  DrawCompare(hdphiBJTdata,hdphiBJTall,"#Delta#phi");
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

  plotdiffmax = 0.15;
  //important! but not used in the result plots
  // DrawCompare(hdtxJASSubEars, hmcxJASSubEars);
  DrawCompare(hdtxJASSubEars, hmcxJASsigBB);


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
  

  // TFile *f = new TFile(Form("xJdphi_bin_%d_%d.root",binMin,binMax),"recreate");
  // hdtINCxJASSubEars->Write(Form("xJ_data_inc_%d_%d",binMin/2,binMax/2));
  // hmcINCxJASSubEars->Write(Form("xJ_mc_inc_%d_%d",binMin/2,binMax/2));
  // hdphiINCdata->Write(Form("dphi_data_inc_%d_%d",binMin/2,binMax/2));
  // hdphiINCall->Write(Form("dphi_mc_inc_%d_%d",binMin/2,binMax/2));
  // f->Close();

}



void tellmetruth(TString name = "", bool applytagg = true, bool tagLJ = true, bool tagSL = true)
{
  name = "results_"+name;
  macro m(name);

  bool applytriggercorr = true;

  tagLeadingJet = tagLJ;

  //tagSL == !mockSL
  dtppjpf = tagSL ? "dtppjpf" : "dXppjpf";
  dtPbbjt = tagSL ? "dtPbbjt" : "dXPbbjt";
  dtPbj60 = tagSL ? "dtPbj60" : "dXPbj60";

 //mcppqcd
 // mcPbbfa
  if (!tagLJ && !tagSL) {    
    mcppqcd = "mXppqcd";
    mcPbqcd = "mXPbqcd";
    mcPbbfa = "mXPbqcd";
  }


  //in case LJ is not tagged - use inc jet ntuple
  if (!tagLJ && tagSL)  dtPbbjt = "dtPbj60";
  if (!tagLJ && !tagSL) dtPbbjt = "dXPbj60";

  applyTriggerCorr = applytriggercorr;
  applyCorrection = applytagg;

//"tellmetruth0503_triggercorrtaggcorr_correctsub"


  loadTagEffCorrections();
  loadTrigEffCorrections();

  findtruthpp(1.0); //fraction of data to process
  findtruthPbPb(60,200);
  findtruthPbPb(20,60); 
  findtruthPbPb(0,20);
  
  //findtruthPbPb(0,200);
  

  cout<<"Bin \t\tInc.Data   \tInc.MC    \tSL.Data   \tSL.MC"<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    cout<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t"<<xjdtincmean[i]<<"\t"<<xjdtincmeanerror[i]<<
                                          "\t"<<xjmcincmean[i]<<"\t"<<xjmcincmeanerror[i]<<
                                          "\t"<<xjdtbjtmean[i]<<"\t"<<xjdtbjtmeanerror[i]<<
                                          "\t"<<xjmcbjtmean[i]<<"\t"<<xjmcbjtmeanerror[i]<<endl;


  std::ofstream ofs (plotfoldername+"/results.txt", std::ofstream::out);

  ofs<<"Bin \t\tInc.Data   \tInc.MC    \tSL.Data   \tSL.MC"<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    ofs<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t"<<xjdtincmean[i]<<"\t"<<xjdtincmeanerror[i]<<
                                          "\t\t"<<xjmcincmean[i]<<"\t"<<xjmcincmeanerror[i]<<
                                          "\t\t"<<xjdtbjtmean[i]<<"\t"<<xjdtbjtmeanerror[i]<<
                                          "\t\t"<<xjmcbjtmean[i]<<"\t"<<xjmcbjtmeanerror[i]<<endl;
  for (unsigned i=0;i<lbinmin.size();i++)
    ofs<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t";
  ofs<<endl;
  for (unsigned i=0;i<xjdtbjtmean.size();i++)
    ofs<<setprecision(3)<<xjdtbjtmean[i]<<"\t";
  ofs<<endl;
    
  ofs.close();


  map<TString, float> res;
  for (unsigned i=0;i<lbinmin.size();i++) {
    const char* end = Form("%d%d",(int)lbinmin[i],(int)lbinmax[i]);
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

  //moneyplot("moneyplot"+suffix);
}
