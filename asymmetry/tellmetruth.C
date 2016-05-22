#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "TRandom.h"
#include "../corrections/tageffcorrections.h"
#include "../corrections/eclipseclosure.C"
#include "moneyplot.C"

// const float pt1cut = 100;
// const float pt2cut = 40;

// const float PI23 = 3.142*2/3;
// const float PI13 = 3.142*1/3;



//vector<float> processWeights(4);

//for reporting in the end
vector<float> lbinmin, lbinmax;
vector<float> xjdtincmean, xjdtincmeanerror;
vector<float> xjdtbjtmean, xjdtbjtmeanerror;

vector<float> xjmcincmean, xjmcincmeanerror;
vector<float> xjmcincsig2mean, xjmcincsig2meanerror;
vector<float> xjmcbjtSBmean,xjmcbjtSBmeanerror;

vector<float> xjmcbjtmean, xjmcbjtmeanerror;

vector<float> xjmcbjtsyslo,xjmcbjtsyshi;

TF1 *ftrigEta, *ftrigCent, *ftrigPt;


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

  return centTagEffCorr;
}


void findtruthpp(float datafraction = 1.)
{
  TFile *fdtpp = new TFile(config.getFileName_djt("dtppjpf"));// "/data_CMS/cms/lisniak/bjet2015/dtppjpfak4PF_djt.root");

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

    if (_["jtpt1"]>pt1cut && _[discr_csvV1_1]>0.9 && _[jtptSL]>pt2cut)
      hdphippBJTdata->Fill(_[dphiSL1],wb);   

    if (_["jtpt1"]>pt1cut && _[discr_csvV1_1]>0.9 && _[jtptSL]>pt2cut && _[dphiSL1]>PI23)
      hdtppxJAS->Fill(_[jtptSL]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && _[discr_csvV1_1]>0.9 && _["jtpt2"]>pt2cut && _[discr_csvV1_2]>0.9 && _["dphi21"]>PI23)
      hdt12ppxJAS->Fill(_["jtpt2"]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && _[discr_csvV1_1]>0.9 && _["jtpt2"]>pt2cut && _[discr_csvV1_2]>0.9)
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

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut)
      hdphippBJTmc->Fill(m[dphiSL1],wb);
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut && m[dphiSL1]>PI23)
      hmcppxJAS->Fill(m[jtptSL]/m["jtpt1"],wb);

    if (m["jtpt1"]>pt1cut && m[discr_csvV1_1]>0.9 && m["jtpt2"]>pt2cut && m[discr_csvV1_2]>0.9)
      hdphippBJT12mc->Fill(m["dphi21"],wb);
    if (m["jtpt1"]>pt1cut && m[discr_csvV1_1]>0.9 && m["jtpt2"]>pt2cut && m[discr_csvV1_2]>0.9 && m["dphi21"]>PI23)
      hmc12ppxJAS->Fill(m["jtpt2"]/m["jtpt1"],wb);


    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcppxJASsigBB->Fill(m["jtptSB"]/m["jtpt1"],wSB);

  });

  TFile *fmcppqcd = new TFile(config.getFileName_djt("mcppqcd"));//("/data_CMS/cms/lisniak/bjet2015/mcppqcdak4PF_djt.root");


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
  aktstring += "R=0.4 |#eta|<2.0";
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
  xjmcbjtmean.push_back(       hmcppxJAS->GetMean());
  xjmcbjtmeanerror.push_back(  hmcppxJAS->GetMeanError());
  xjmcbjtSBmean.push_back(         hmcppxJASsigBB->GetMean());//hmcxJASsigBB->GetMean());
  xjmcbjtSBmeanerror.push_back(    hmcppxJASsigBB->GetMeanError());//hmcxJASsigBB->GetMeanError());

  cout<<"Is it 0 ? "<<hmcppqcdxJASsignal2->GetMean()<<endl;

  xjmcincsig2mean.push_back(        hmcppqcdxJASsignal2->GetMean());
  xjmcincsig2meanerror.push_back(  hmcppqcdxJASsignal2->GetMeanError());

  xjmcbjtsyslo.push_back(0.001);
  xjmcbjtsyshi.push_back(0.001);

  TFile *f = new TFile("xJdphi_pp.root","recreate");
  hdtINCppxJAS->Write("xJ_data_inc_pp");
  hmcppqcdxJAS->Write("xJ_mc_inc_pp");
  f->Close();

}

void findtruthPbPb(int binMin, int binMax)
{
  TString dtfname = config.getFileName_djt("dtPbbjt");
  TString mcfname = config.getFileName_djt("mcPbbfa");

  TFile *fdt = new TFile(dtfname);
  TFile *fmc = new TFile(mcfname);
  TFile *fdtinc = new TFile(config.getFileName_djt("dtPbj60"));
  TFile *fmcinc = new TFile(config.getFileName_djt("mcPbqcd"));



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

  // auto hdtxJASSub = geth("hdtxJASSub","SL Data SB subtracted;x_{J}");
  // auto hmcxJASSub = geth("hmcxJASSub","SL MC SB subtracted;x_{J}");
  auto hdtxJASSubEars = geth("hdtxJASSubEars","Data b-jets;x_{J}");
  auto hmcxJASSubEars = geth("hmcxJASSubEars","MC b-jets;x_{J}");

  auto hmcxJASsig = geth("hmcxJASsig","SL sig MC;x_{J}");
  auto hmcxJASsigsyslo = geth("hmcxJASsigsyslo","SL sig MC sys lo;x_{J}");
  auto hmcxJASsigsyshi = geth("hmcxJASsigsyshi","SL sig MC sys hi;x_{J}");
  auto hmcxJASsigBB = geth("hmcxJASsigBB","MC b-jets;x_{J}");


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

  Fill(fdt,{"weight","vz","bin","jtpt1",discr_csvV1_1,jtptSL,dphiSL1,"jteta1",jtetaSL,"numTagged"},[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL],m["bin"]);
    float trigcorr = getTrigcorrection(m["jtpt1"],m["jteta1"],m["bin"]); //WRONG! use trigger pt
    float wb = w*corr*trigcorr;
    float dphi = m[dphiSL1];
    float deta = abs(m["jteta1"]-m[jtetaSL]);
    hdtvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut) {

      hdtbin->Fill(m["bin"],wb);
      hdphiBJTdata->Fill(m[dphiSL1],wb);

      if (dphi>PI23)
        hdtxJAS->Fill(m[jtptSL]/m["jtpt1"],wb);
      if (dphi<PI13)
        hdtxJNS->Fill(m[jtptSL]/m["jtpt1"],wb);

      if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
        hdtxJEars->Fill(m[jtptSL]/m["jtpt1"],wb);

    }

  });

  vector<float> cbins = {0.,20.,60.,200.};
  vector<TString> binnames = {"0-10%", "10-30%", "30-100%"};
  int Nb = 3;
  vector<TH1F *>hsig(Nb);
  vector<TH1F *>hasd(Nb);
  vector<TH1F *>hbkg(Nb);
  vector<TH1F *>hsub(Nb);
  vector<TH1F *>hhyj(Nb);
  vector<TH1F *>hshj(Nb);
  // vector<TH1F *>hsigsyslo(Nb);
  // vector<TH1F *>hsigsyshi(Nb);

  vector<TH1F *>hptNSsig(Nb);
  vector<TH1F *>hptNSbkg(Nb);
  vector<TH1F *>hptNS(Nb);

  for (int i=0;i<Nb;i++) {
    seth(10,0,1);//cbins);//10,0,1);
    hsig[i] = geth(Form("hsig%d",i),Form("Signal away-side %s;x_{J}",binnames[i].Data())) ;
    hasd[i] = geth(Form("hasd%d",i),Form("Measured away-side %s;x_{J}",binnames[i].Data()));
    hbkg[i] = geth(Form("hbkg%d",i),Form("Near-side %s;x_{J}",binnames[i].Data()));
    hhyj[i] = geth(Form("hhyj%d",i),Form("dphi<1/3pi hydjet %s;x_{J}",binnames[i].Data()));
    hsub[i] = geth(Form("hsub%d",i),Form("Subtracted NS %s;x_{J}",binnames[i].Data()));
    hshj[i] = geth(Form("hshj%d",i),Form("Subtracted Hydjet %s;x_{J}",binnames[i].Data()));

    // hsigsyslo[i] = geth(Form("hsigsyslo%d",i),Form("Signal sys low %s;x_{J}",binnames[i].Data())) ;
    // hsigsyshi[i] = geth(Form("hsigsyshi%d",i),Form("Signal sys high %s;x_{J}",binnames[i].Data())) ;


    seth(10,40,100);
    hptNSsig[i] = geth(Form("hptNSsig%d",i),Form("Near-side signal %s;p_{T} [GeV]",binnames[i].Data()));
    hptNS[i] = geth(Form("hptNS%d",i),Form("Near-side total %s;p_{T} [GeV]",binnames[i].Data()));
    hptNSbkg[i] = geth(Form("hptNSbkg%d",i),Form("Near-side bkg %s;p_{T} [GeV]",binnames[i].Data()));
  }

  int Nbinc = 10;
  vector<TH1F *>hincsig(Nbinc);
  vector<TH1F *>hincasd(Nbinc);
  vector<TH1F *>hincbkg(Nbinc);
  vector<TH1F *>hincsub(Nbinc);
  for (int i=0;i<Nbinc;i++) {
    hincsig[i] = geth(Form("hincsig%d",i));
    hincasd[i] = geth(Form("hincasd%d",i));
    hincbkg[i] = geth(Form("hincbkg%d",i));
    hincsub[i] = geth(Form("hincsub%d",i));
  }



  Fill(fmc,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1",discr_csvV1_1,jtptSL,dphiSL1,"jteta1",jtetaSL,subidSL,refptSL,pairCodeSL1,
            "jtptSignal2",discr_csvV1_Signal2,"Signal2ord","SLord","dphiSignal21",refparton_flavorForBSL,"refparton_flavorForB1","numTagged","jtptSB","dphiSB1"},[&] (dict &m) {
    if (m["numTagged"]>6) return;
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["pthat"]<pthatcut) return; /////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //at least one of the two jets must be a b-jet
    if (abs(m["refparton_flavorForB1"])!=5 && abs(m[refparton_flavorForBSL])!=5) return;

    //float w0 = m["weight"];


    float dphi = m[dphiSL1];
    float deta = abs(m["jteta1"]-m[jtetaSL]);

    float w0=weight1SLPbPb(m);//w0*processWeights[(int)m["bProdCode"]];

    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m[jtptSL],m[jtetaSL],m["bin"]);
    float w = w0*corr;

    float wSB = m["weight"]*processweight((int)m["bProdCode"]);

    //float wbkg = w0*corr;//*(m["pthat"]>80);

    hmcvz->Fill(m["vz"],w);


    //do that only in 0-100% pass
    if (binMin==0 && binMax==200) {
      //int i=((int)m["bin"])/20;
      int i=0;
      int bin = m["bin"];
      if (bin<20) i=0;
      if (bin>=20 && bin<60) i=1;
      if (bin>=60) i=2;


//&& abs(m["refparton_flavorForB1"])==5
//&& m[discr_csvV1_1]>0.9

//signal: subidSL==0 on the away side
//background: subidSL!=0 on the away side

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m[jtptSL]>pt2cut && m[dphiSL1]>PI23)
        hasd[i]->Fill(m[jtptSL]/m["jtpt1"], w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);//wbkg); ///////WRONG POTENTIALLY!

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m[jtptSL]>pt2cut && m[dphiSL1]<PI13 && !IsSignal(m))
        hhyj[i]->Fill(m[jtptSL]/m["jtpt1"], w);//wbkg);

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m[jtptSL]>pt2cut && m[dphiSL1]>PI23 && IsSignal(m)) {//&& m[pairCodeSL1]==0) {
          hsig[i]->Fill(m[jtptSL]/m["jtpt1"], w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);
          // float pt2lo = gRandom->Uniform() < 0.1 ? pt2cut     : m[jtptSL];
          // float pt2hi = gRandom->Uniform() > 0.9 ? m["jtpt1"] : m[jtptSL];
          // hsigsyslo[i]->Fill(pt2lo/m["jtpt1"], w);
          // hsigsyshi[i]->Fill(pt2hi/m["jtpt1"], w);
      }
      

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m[jtptSL]>pt2cut
//           && m[dphiSL1]<PI13) {
           && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))) {
        hbkg[i]->Fill(m[jtptSL]/m["jtpt1"],w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);//wbkg);?????????????
        hptNS[i]->Fill(m[jtptSL], w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);
        if (IsSignal(m))
          hptNSsig[i]->Fill(m[jtptSL], w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);
        else
          hptNSbkg[i]->Fill(m[jtptSL], w);//abs(m[refparton_flavorForBSL])==5  ?  wb : wbkg);
      }

    }



//NOT CORRECTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m[discr_csvV1_1]>0.9 && m["jtptSignal2"]>pt2cut && m[discr_csvV1_Signal2]>0.9) { //&& m["dphiSignal21"]>PI23 
      hbinSignal->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==2)
        hbinSignalFound12->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==m[SLord])
        hbinSignalFoundSL->Fill(m["bin"]/2,w);
    }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut && m[dphiSL1]>PI23) {
      hbinSL->Fill(m["bin"],w);//wb);
      if (abs(m[refparton_flavorForBSL])==5)
        hbinSLisB->Fill(m["bin"],w);//wb);
    }


    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcxJASsigBB->Fill(m["jtptSB"]/m["jtpt1"],wSB);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut) {
      hmcbin->Fill(m["bin"],w);//wb);
      hdphiBJTall->Fill(m[dphiSL1],w);//wb);

      if (m[pairCodeSL1]==0)
        hdphiBJTsig->Fill(m[dphiSL1],w);//wb);

      if (m[dphiSL1]>PI23) {
        hmcxJAS->Fill(m[jtptSL]/m["jtpt1"], w);//m[pairCodeSL1]==0 ? wb : wbkg);

        //signal xJ
        if (m[pairCodeSL1]==0) {
          hmcxJASsig->Fill(m[jtptSL]/m["jtpt1"],w);//wb);
          float pt2lo = gRandom->Uniform() < 0.05 ? pt2cut     : m[jtptSL];
          float pt2hi = gRandom->Uniform() > 0.95 ? m["jtpt1"] : m[jtptSL];
          hmcxJASsigsyslo->Fill(pt2lo/m["jtpt1"], w);
          hmcxJASsigsyshi->Fill(pt2hi/m["jtpt1"], w);
        }

        //signal-like pairCode (for MC purity)
        if (IsSignal(m))//m[subidSL]==0)
          hPairCodeBFA->Fill(m[pairCodeSL1],w0);//wb);////////////////////////////////there were w0 here
      }

      if (m[dphiSL1]<PI13)
        hmcxJNS->Fill(m[jtptSL]/m["jtpt1"],w);//wb);


    }

// // here I subtract background with TRUE leading jet!!!

// float LJeff = 0.35;

if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m[discr_csvV1_1]>0.9  && m[jtptSL]>pt2cut //&& abs(m["refparton_flavorForB1"])==5
        && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1)))
        hmcxJEars->Fill(m[jtptSL]/m["jtpt1"],w);//LFeff*w//LJeff*(m[pairCodeSL1]==0 ? wb : wbkg));

  });



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


   //check closure that only in 0-100% pass
   if (binMin==0 && binMax==200) {
    int i=((int)m["bin"])/20;
    //int i=0;
    //int bin = m["bin"];
    //if (bin<20) i=0;
    //if (bin>=20 && bin<60) i=1;
    //if (bin>=60) i=2;


    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut && m["dphi21"]>PI23)
      hincasd[i]->Fill(m["jtpt2"]/m["jtpt1"],w);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut && m["dphi21"]>PI23 && m["subid2"]==0)
      hincsig[i]->Fill(m["jtpt2"]/m["jtpt1"],w);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut
         && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1)))
      hincbkg[i]->Fill(m["jtpt2"]/m["jtpt1"],w);

   }




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

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m[discr_csvV1_1]>0.9 && m[jtptSL]>pt2cut && m[dphiSL1]>PI23) {
      if (IsSignal(m) && m["numTagged"]<=6)//m[pairCodeSL1]<4 && 
        hPairCodeQCD->Fill(m[pairCodeSL1],w);
    }

  });


  hPairCode->SetBinContent(1,hPairCodeBFA->GetBinContent(1));
  hPairCode->SetBinContent(2,hPairCodeQCD->GetBinContent(2));
  hPairCode->SetBinContent(3,hPairCodeBFA->GetBinContent(3));
  hPairCode->SetBinContent(4,hPairCodeBFA->GetBinContent(4));
  hPairCode->SetBinContent(5,hPairCodeQCD->GetBinContent(5));

//  hPairCode->SetBinContent(5,hPairCodeBFA->GetBinContent(5)*30); //cheating, interesting idea

  // Normalize({hdtvz,hmcvz,hdtbin,hmcbin});

  // Draw({hdtvz,hmcvz});
  // plotylog = true;
  // Draw({hdtbin,hmcbin});

  //hdtxJASSub->Add(hdtxJAS,hdtxJNS,1,-1);
  //hmcxJASSub->Add(hmcxJAS,hmcxJNS,1,-1);


  float syslo = hmcxJASsig->GetMean() - hmcxJASsigsyslo->GetMean();
  float syshi = hmcxJASsigsyshi->GetMean() - hmcxJASsig->GetMean();

//  cout<<"sys lo/hi "<<binMin<<" "<<binMax<<syslo<<" : "<<syshi<<endl;

  //TODO: fix bin handling
  float coef = bkgfractionInNearSide[getbinindex((binMin+binMax)/2)];
  hdtxJASSubEars->Add(hdtxJAS,hdtxJEars,1,-1*coef);
  hmcxJASSubEars->Add(hmcxJAS,hmcxJEars,1,-1*coef);

  cout<<"   Fraction of data left after subtraction in bin "<<binMin/2<<" - "<<binMax/2<<" = "<<hdtxJASSubEars->Integral()/hdtxJAS->Integral()<<endl;

  hdtINCxJASSubEars->Add(hdtINCxJAS,hdtINCxJEars,1,-1);
  hmcINCxJASSubEars->Add(hmcINCxJAS,hmcINCxJEars,1,-1);





  for (int i=0;i<Nb;i++) {
    hsub[i]->Add(hasd[i],hbkg[i],1,-1*bkgfractionInNearSide[i]);
    hshj[i]->Add(hasd[i],hhyj[i],1,-1);
  }
  for (int i=0;i<Nbinc;i++) 
    hincsub[i]->Add(hincasd[i],hincbkg[i],1,-1);

  seth(cbins);//Nb,0,100);
  auto hcentrSubSIG = geth("hcentrSubSIG","Signal;bin;<x_{J}>");
  auto hcentrSubASD = geth("hcentrSubASD","Unsubtracted;bin;<x_{J}>");
  auto hcentrSubBKG = geth("hcentrSubBKG","Background;bin;<x_{J}>");
  auto hcentrSubCLS = geth("hcentrSubCLS","Subtracted;bin;<x_{J}>");
  auto hcentrSubHJS = geth("hcentrSubHJS","Subtracted Hydjet;bin;<x_{J}>");


  seth(Nbinc,0,100);
  auto hcentrSubSIGinc = geth("hcentrSubSIGinc","INC Signal;bin;<x_{J}>");
  auto hcentrSubASDinc = geth("hcentrSubASDinc","INC Unsubtracted;bin;<x_{J}>");
  auto hcentrSubBKGinc = geth("hcentrSubBKGinc","INC Background;bin;<x_{J}>");
  auto hcentrSubCLSinc = geth("hcentrSubCLSinc","INC Subtracted;bin;<x_{J}>"); 

  for (int i=0;i<Nb;i++) {
    hcentrSubSIG->SetBinContent(i+1,hsig[i]->GetMean());hcentrSubSIG->SetBinError(i+1,hsig[i]->GetMeanError());
    hcentrSubASD->SetBinContent(i+1,hasd[i]->GetMean());hcentrSubASD->SetBinError(i+1,hasd[i]->GetMeanError());
    hcentrSubBKG->SetBinContent(i+1,hbkg[i]->GetMean());hcentrSubBKG->SetBinError(i+1,hbkg[i]->GetMeanError());

    hcentrSubCLS->SetBinContent(i+1,hsub[i]->GetMean());hcentrSubCLS->SetBinError(i+1,hsub[i]->GetMeanError());
    hcentrSubHJS->SetBinContent(i+1,hshj[i]->GetMean());hcentrSubHJS->SetBinError(i+1,hshj[i]->GetMeanError());
  }
  Print(hcentrSubBKG);
  for (int i=0;i<Nbinc;i++) {
    hcentrSubSIGinc->SetBinContent(i+1,hincsig[i]->GetMean());hcentrSubSIGinc->SetBinError(i+1,hincsig[i]->GetMeanError());
    hcentrSubASDinc->SetBinContent(i+1,hincasd[i]->GetMean());hcentrSubASDinc->SetBinError(i+1,hincasd[i]->GetMeanError());
    hcentrSubBKGinc->SetBinContent(i+1,hincbkg[i]->GetMean());hcentrSubBKGinc->SetBinError(i+1,hincbkg[i]->GetMeanError());
    hcentrSubCLSinc->SetBinContent(i+1,hincsub[i]->GetMean());hcentrSubCLSinc->SetBinError(i+1,hincsub[i]->GetMeanError());
  }



//hcentrSubASD->Add(hcentrSubSIG,-1);
//hcentrSubBKG->Add(hcentrSubSIG,-1);
//hcentrSubCLS->Add(hcentrSubSIG,-1);
//hcentrSubHJS->Add(hcentrSubSIG,-1);
//hcentrSubSIG->Add(hcentrSubSIG,-1);





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
  aktstring = "anti-k_{T} Pu R=0.4 |#eta|<2.0";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2b}>%d GeV", (int)pt1cut, (int)pt2cut);
  
  plotytitle = "Event fractions";

  plottextposx = 10;
  plotlegendpos = TopRight;
  plotputmean = false;
  plotymin = 0.;

  // for (int i=0;i<Nb;i++) {
  //         plotymax = 1E-8;
  //   Draw({hbkg[i],hhyj[i]});
  //         plotymax = 1E-7;
  //   Draw({hsig[i],hasd[i],hsub[i],hshj[i]});
  // }

for (int i=0;i<Nb;i++)
  cout<<"Amount of background in the subtraction "<<binnames[i]<<" : "<<
      setprecision(2)<<(1-hptNSsig[i]->Integral(0,hptNSsig[i]->GetNbinsX()+1)/(float)hptNS[i]->Integral(0,hptNS[i]->GetNbinsX()+1))*100<<" % "<<endl;;



 //Normalize({hptNSsig[0],hptNSsig[1],hptNSsig[2],hptNS[0],hptNS[1],hptNS[2]});
//   Draw({hptNSsig[0],hptNSsig[1],hptNSsig[2]});
//   Draw({hptNS[0],hptNS[1],hptNS[2]});


for (int i=0;i<Nb;i++) {
  //subtract signal from all to get bkg only
  float norm = hptNSbkg[i]->Integral()/hptNSsig[i]->Integral();
  hptNSsig[i]->Scale(1/hptNSsig[i]->Integral());
  hptNSbkg[i]->Scale(norm/hptNSbkg[i]->Integral());
  
  plotymax = 1.;//2.5E-8;
  // Draw({hptNSsig[i],hptNSbkg[i]});
}

  plotymin = 0.5;//0.4;
  plotymax = 0.8;//0.8;
  plotlegendpos = BottomRight;
  plottextposx = 0.5;
  plottextposy = 0.79;

  plotputmean = false;
  // Draw({hcentrSubSIG, hcentrSubASD, hcentrSubCLS,hcentrSubHJS});



  // Draw({hcentrSubSIGinc, hcentrSubASDinc, hcentrSubCLSinc});



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




//oops
  
  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcxJASsigBB,hmcINCxJASSubEars,hmcINCxJASsignal2});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsigBB,hmcxJASsig,hmcxJASsigsyshi,hmcxJASsigsyslo});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsignal2,hmcINCxJASsig});





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
  // Draw({hPairCode});


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
  xjmcbjtmean.push_back(       hmcxJASSubEars->GetMean());
  xjmcbjtmeanerror.push_back(  hmcxJASSubEars->GetMeanError());
  xjmcbjtSBmean.push_back(         hmcxJASsigBB->GetMean());
  xjmcbjtSBmeanerror.push_back(    hmcxJASsigBB->GetMeanError());
  xjmcincsig2mean.push_back(       hmcINCxJASsignal2->GetMean());
  xjmcincsig2meanerror.push_back(  hmcINCxJASsignal2->GetMeanError());
  

  cout<<"sys lo/hi "<<binMin<<" "<<binMax<<" "<<syslo<<" : "<<syshi<<endl;

  xjmcbjtsyslo.push_back(syslo);
  xjmcbjtsyshi.push_back(syshi);

  TFile *f = new TFile(Form("xJdphi_bin_%d_%d.root",binMin,binMax),"recreate");
  hdtINCxJASSubEars->Write(Form("xJ_data_inc_%d_%d",binMin/2,binMax/2));
  hmcINCxJASSubEars->Write(Form("xJ_mc_inc_%d_%d",binMin/2,binMax/2));
  hdphiINCdata->Write(Form("dphi_data_inc_%d_%d",binMin/2,binMax/2));
  hdphiINCall->Write(Form("dphi_mc_inc_%d_%d",binMin/2,binMax/2));
  f->Close();

}



void tellmetruth(TString name = "")
{
  name = "results_"+name;
  macro m(name);

  int subtractcode = 2;
  bool applytriggercorr = true;
  bool applytagg = true;




  if (subtractcode==0) bkgfractionInNearSide = {0,0,0};
  if (subtractcode==1) bkgfractionInNearSide = {1,1,1};
//subtractcode = 2 or something - use whatever is true

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

    res[Form("xj_mc_bjt_mean%s",end)]      = xjmcbjtmean[i];
    res[Form("xj_mc_bjt_meanerror%s",end)] = xjmcbjtmeanerror[i];

    res[Form("xj_mc_bjtSB_mean%s",end)]      = xjmcbjtSBmean[i];
    res[Form("xj_mc_bjtSB_meanerror%s",end)] = xjmcbjtSBmeanerror[i];


    cout<<Form("xjmcbjtsyslo%s",end)<<xjmcbjtsyslo[i]<<endl;
    res[Form("xjmcbjtsyslo%s",end)]      = xjmcbjtsyslo[i];
    res[Form("xjmcbjtsyshi%s",end)] = xjmcbjtsyshi[i];

  }

  WriteToFile(plotfoldername+"/results.root",res);

  //moneyplot("moneyplot"+suffix);
}
