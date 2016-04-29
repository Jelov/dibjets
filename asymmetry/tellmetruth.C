#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

// const float pt1cut = 100;
// const float pt2cut = 40;

// const float PI23 = 3.142*2/3;
// const float PI13 = 3.142*1/3;

bool applyCorrection = false;

//vector<float> processWeights(4);

//for reporting in the end
vector<float> lbinmin, lbinmax;
vector<float> xjdtincmean, xjdtincmeanerror;
vector<float> xjdtbjtmean, xjdtbjtmeanerror;

vector<float> xjmcincmean, xjmcincmeanerror;
vector<float> xjmcbjtmean, xjmcbjtmeanerror;

TF1 *fppEta, *fppPt, *fPbPbCent, *fPbPbEta, *fPbPbPt;


void loadTagEffCorrections()
{
   TFile *fppFits = new TFile("/data_CMS/cms/mnguyen/bJet2015/tagEffCorr/offlineTagEff_PP.root");
   
   fppEta = (TF1*) fppFits->Get("fitEta");
   fppPt = (TF1*) fppFits->Get("fitPt");


  TFile *fPbPbFits = new TFile("/data_CMS/cms/mnguyen/bJet2015/tagEffCorr/offlineTagEff_PbPb.root");

  fPbPbCent = (TF1*) fPbPbFits->Get("fitCent");
  fPbPbEta = (TF1*) fPbPbFits->Get("fitEta");
  fPbPbPt = (TF1*) fPbPbFits->Get("fitPt");
}

float getppcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta)
{

  if (!applyCorrection) return 1;

  double leadingPtTagEffCorr = (leadingjetpt > 200) ? fppPt->Eval(200):fppPt->Eval(leadingjetpt);
  double leadingTagEffCorr = 1./fppEta->Eval(leadingjeteta)/leadingPtTagEffCorr;


  double parterPtTagEffCorr = (partnerjetpt > 200) ? fppPt->Eval(200):fppPt->Eval(partnerjetpt);
  double partnerTagEffCorr = 1./fppEta->Eval(partnerjeteta)/parterPtTagEffCorr;
     
  // numerical factor will just give an average correction of 1, can omit this
  double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr/8.98489;
  return combinedTagEffCorr;
}

float getPbPbcorrection(float leadingjetpt, float leadingjeteta, float partnerjetpt, float partnerjeteta, float bin)
{

  if (!applyCorrection) return 1;

  double centTagEffCorr = 1./fPbPbCent->Eval(bin/2.)/fPbPbCent->Eval(bin/2.);

  double leadingPtTagEff = (leadingjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(leadingjetpt);
  double leadingTagEffCorr = 1./fPbPbEta->Eval(leadingjeteta)/leadingPtTagEff;


  double parterPtTagEff = (partnerjetpt > 200) ? fPbPbPt->Eval(200):fPbPbPt->Eval(partnerjetpt);
  double partnerTagEffCorr = 1./fPbPbEta->Eval(partnerjeteta)/parterPtTagEff;

  double combinedTagEffCorr =leadingTagEffCorr*partnerTagEffCorr*centTagEffCorr/488.652;
  return combinedTagEffCorr;
}

void findtruthpp()
{
  TFile *fdtpp = new TFile(config.getFileName_djt("dtppjpf"));// "/data_CMS/cms/lisniak/bjet2015/dtppjpfak4PF_djt.root");

  buildh(10,0,1);
  auto hdtppxJAS = geth("hdtppxJAS","Data b-jets SL;x_{J}");
  auto hdtINCppxJAS = geth("hdtINCppxJAS","Data Inclusive;x_{J}");
  auto hdt12ppxJAS = geth("hdt12ppxJAS","Data b-jets 12;x_{J}");

  auto hmcppxJAS = geth("hmcppxJAS","MC b-jets SL;x_{J}");
  auto hmcppqcdxJAS = geth("hmcppqcdxJAS","MC Inclusive;x_{J}");
  auto hmc12ppxJAS = geth("hmc12ppxJAS","MC b-jets 12;x_{J}");

  Fill(fdtpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jtpt2","dphi21","discr_csvV1_2","jteta1","jtetaSL"},[&] (dict _) {
    float w = _["weight"];
    float corr = getppcorrection(_["jtpt1"],_["jtpt2"],_["jtptSL"],_["jtetaSL"]);
    float wb = w*corr;

    if (_["jtpt1"]>pt1cut && _["discr_csvV1_1"]>0.9 && _["jtptSL"]>pt2cut && _["dphiSL1"]>PI23)
      hdtppxJAS->Fill(_["jtptSL"]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && _["discr_csvV1_1"]>0.9 && _["jtpt2"]>pt2cut && _["discr_csvV1_2"]>0.9 && _["dphi21"]>PI23)
      hdt12ppxJAS->Fill(_["jtpt2"]/_["jtpt1"],wb);

    if (_["jtpt1"]>pt1cut && _["jtpt2"]>pt2cut && _["dphi21"]>PI23)
      hdtINCppxJAS->Fill(_["jtpt2"]/_["jtpt1"],w);

    },0.2);

  TFile *fmcpp = new TFile(config.getFileName_djt("mcppbfa"));//("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");
  Fill(fmcpp,{"pthat","weight","jtpt1","refpt1","discr_csvV1_1","jtptSL","dphiSL1","bProdCode","jtpt2", "discr_csvV1_2", "dphi21","jtetaSL"},[&] (dict &m) {
    if (m["pthat"]<65) return;
    float w = m["weight"]*processWeights[(int)m["bProdCode"]];
    float corr = getppcorrection(m["jtpt1"],m["jtpt2"],m["jtptSL"],m["jtetaSL"]);
    float wb = w*corr;

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>PI23)
      hmcppxJAS->Fill(m["jtptSL"]/m["jtpt1"],wb);
    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]>PI23)
      hmc12ppxJAS->Fill(m["jtpt2"]/m["jtpt1"],wb);

  });

  TFile *fmcppqcd = new TFile(config.getFileName_djt("mcppqcd"));//("/data_CMS/cms/lisniak/bjet2015/mcppqcdak4PF_djt.root");


  Fill(fmcppqcd,{"pthat","weight","jtpt1","refpt1","jtpt2","dphi21"},[&] (dict &m) {
    if (m["pthat"]<65) return;
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtpt2"]>pt2cut && m["dphi21"]>PI23) {
        hmcppqcdxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

  });

  SetData({hdtppxJAS,hdtINCppxJAS,hdt12ppxJAS});
  SetMC({hmcppxJAS,hmcppqcdxJAS,hmc12ppxJAS});
  SetInc({hdtINCppxJAS,hmcppqcdxJAS});
  SetB({hdtppxJAS,hmcppxJAS,hdt12ppxJAS,hmc12ppxJAS});

  NormalizeAllHists();
  plotputmean = true;
  plotytitle = "Event fractions";
  plotdivide = false;
  aktstring += "R=0.4 |#eta|<2.0";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";
  plotymax = 0.4;

  DrawCompare(hdtppxJAS,hmcppxJAS);
  DrawCompare(hdtINCppxJAS,hmcppqcdxJAS);
  DrawCompare(hdt12ppxJAS,hmc12ppxJAS); 


  lbinmin.push_back(-1);lbinmax.push_back(-1);
  xjdtincmean.push_back(       hdtINCppxJAS->GetMean());
  xjdtincmeanerror.push_back(  hdtINCppxJAS->GetMeanError());
  xjmcincmean.push_back(       hmcppqcdxJAS->GetMean());
  xjmcincmeanerror.push_back(  hmcppqcdxJAS->GetMeanError());
  xjdtbjtmean.push_back(       hdtppxJAS->GetMean());
  xjdtbjtmeanerror.push_back(  hdtppxJAS->GetMeanError());
  xjmcbjtmean.push_back(       hmcppxJAS->GetMean());
  xjmcbjtmeanerror.push_back(  hmcppxJAS->GetMeanError());
}

//bool IsSignal(dict d) { return d["subidSL"]==0 && d["refptSL"]>20;}


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
  buildh(40,-15,15);
  auto hdtvz = geth("hdtvz","Data;vz [cm]");
  auto hmcvz = geth("hmcvz","MC;vz [cm]");

  buildh(100,0,200);
  auto hdtbin = geth("hdtbin","Data;bin");
  auto hmcbin = geth("hmcbin","MC;bin"); 

  //xJ
  buildh(10,0,1);
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


  //inclusive jets
  auto hdtINCxJAS = geth("hdtINCxJAS","INC Data away-side;x_{J}");
  auto hdtINCxJEars = geth("hdtINCxJEars","INC Data ears;x_{J}");
  auto hmcINCxJAS = geth("hmcINCxJAS","INC MC away-side;x_{J}");
  auto hmcINCxJEars = geth("hmcINCxJEars","INC MC ears;x_{J}");

  auto hmcINCxJASsig = geth("hmcINCxJASsig","INC sig MC;x_{J}");


  auto hdtINCxJASSubEars = geth("hdtINCxJASSubEars","Data Inclusive;x_{J}");
  auto hmcINCxJASSubEars = geth("hmcINCxJASSubEars","MC Inclusive;x_{J}");




  //dphi
  buildh(20,0,3.142);
  auto hdphiINCdata = geth("hdphiINCdata","Data Inclusive;#Delta#phi");
  auto hdphiINCall = geth("hdphiINCall","MC Inclusive;#Delta#phi");
  auto hdphiINCsig = geth("hdphiINCsig","MC Inclusive, signal;#Delta#phi");
  auto hdphiBJTdata = geth("hdphiBJTdata","Data b-jets;#Delta#phi");
  auto hdphiBJTall = geth("hdphiBJTall","MC b-jets;#Delta#phi");
  auto hdphiBJTsig = geth("hdphiBJTsig","MC b-jets, signal;#Delta#phi");



  //pair codes
  buildh(5,0,5);
  auto hPairCodeQCD = geth("hPairCodeQCD");
  auto hPairCodeBFA = geth("hPairCodeBFA");
  auto hPairCode = geth("hPairCode");

  //centrality
  buildh(10,0,100); //don't forget to divide!!!
  auto hbinconfusion12 = geth("hbinconfusion12","12 analysis;bin;Jet Confusion");
  auto hbinconfusionSL = geth("hbinconfusionSL","SL analysis;bin;Jet Confusion");
  auto hbinSignal = geth("hbinSignal");
  auto hbinSignalFound12 = geth("hbinSignalFound12");
  auto hbinSignalFoundSL = geth("hbinSignalFoundSL");

  auto hbinfakerateSL = geth("hbinfakerateSL","SL analysis;bin;Purity");
  auto hbinSL = geth("hbinSL");
  auto hbinSLisB = geth("hbinSLisB");

  Fill(fdt,{"weight","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m["jtptSL"],m["jtetaSL"],m["bin"]);
    float wb = w*corr;
    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);
    hdtvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {

      hdtbin->Fill(m["bin"],wb);
      hdphiBJTdata->Fill(m["dphiSL1"],wb);

      if (dphi>PI23)
        hdtxJAS->Fill(m["jtptSL"]/m["jtpt1"],wb);
      if (dphi<PI13)
        hdtxJNS->Fill(m["jtptSL"]/m["jtpt1"],wb);

      if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
        hdtxJEars->Fill(m["jtptSL"]/m["jtpt1"],wb);

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
  vector<TH1F *>hptNSsig(Nb);
  vector<TH1F *>hptNS(Nb);

  for (int i=0;i<Nb;i++) {
    buildh(10,0,1);//cbins);//10,0,1);
    hsig[i] = geth(Form("hsig%d",i),Form("Signal away-side %s;x_{J}",binnames[i].Data())) ;
    hasd[i] = geth(Form("hasd%d",i),Form("Measured away-side %s;x_{J}",binnames[i].Data()));
    hbkg[i] = geth(Form("hbkg%d",i),Form("Near-side %s;x_{J}",binnames[i].Data()));
    hhyj[i] = geth(Form("hhyj%d",i),Form("dphi<1/3pi hydjet %s;x_{J}",binnames[i].Data()));
    hsub[i] = geth(Form("hsub%d",i),Form("Subtracted NS %s;x_{J}",binnames[i].Data()));
    hshj[i] = geth(Form("hshj%d",i),Form("Subtracted Hydjet %s;x_{J}",binnames[i].Data()));
    buildh(10,40,100);
    hptNSsig[i] = geth(Form("hptNSsig%d",i),Form("Near-side signal %s;p_{T} GeV",binnames[i].Data()));
    hptNS[i] = geth(Form("hptNS%d",i),Form("Near-side %s;p_{T} GeV",binnames[i].Data()));
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



  Fill(fmc,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL","subidSL","refptSL","pairCodeSL1",
            "jtptSignal2","discr_csvV1_Signal2","Signal2ord","SLord","dphiSignal21","refparton_flavorForBSL","refparton_flavorForB1"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["pthat"]<80) return; /////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //at least one of the two jets must be a b-jet
    if (abs(m["refparton_flavorForB1"])!=5 && abs(m["refparton_flavorForBSL"])!=5) return;

    float w0 = m["weight"];


    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);

    float w=w0*processWeights[(int)m["bProdCode"]];

    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m["jtptSL"],m["jtetaSL"],m["bin"]);
    float wb = w*corr;


    float wbkg = w0*corr;//*(m["pthat"]>80);

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
//&& m["discr_csvV1_1"]>0.9

//signal: subidSL==0 on the away side
//background: subidSL!=0 on the away side

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSL"]>pt2cut && m["dphiSL1"]>PI23)
        hasd[i]->Fill(m["jtptSL"]/m["jtpt1"], abs(m["refparton_flavorForBSL"])==5  ?  wb : wbkg);//wbkg); ///////WRONG POTENTIALLY!

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSL"]>pt2cut && m["dphiSL1"]<PI13 && !IsSignal(m))
        hhyj[i]->Fill(m["jtptSL"]/m["jtpt1"],wbkg);

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSL"]>pt2cut && m["dphiSL1"]>PI23 && IsSignal(m)) {//&& m["pairCodeSL1"]==0) {
          hsig[i]->Fill(m["jtptSL"]/m["jtpt1"], abs(m["refparton_flavorForBSL"])==5  ?  wb : wbkg);

      }
      

      if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSL"]>pt2cut
//           && m["dphiSL1"]<PI13) {
           && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))) {
        hbkg[i]->Fill(m["jtptSL"]/m["jtpt1"],abs(m["refparton_flavorForBSL"])==5  ?  wb : wbkg);//wbkg);?????????????
        hptNS[i]->Fill(m["jtptSL"], abs(m["refparton_flavorForBSL"])==5  ?  wb : wbkg);
        if (IsSignal(m))
          hptNSsig[i]->Fill(m["jtptSL"], abs(m["refparton_flavorForBSL"])==5  ?  wb : wbkg);
      }

    }



//NOT CORRECTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["discr_csvV1_1"]>0.9 && m["jtptSignal2"]>pt2cut && m["discr_csvV1_Signal2"]>0.9) { //&& m["dphiSignal21"]>PI23 
      hbinSignal->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==2)
        hbinSignalFound12->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==m["SLord"])
        hbinSignalFoundSL->Fill(m["bin"]/2,w);
    }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>PI23) {
      hbinSL->Fill(m["bin"],wb);
      if (abs(m["refparton_flavorForBSL"])==5)
        hbinSLisB->Fill(m["bin"],wb);
    }



    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {
      hmcbin->Fill(m["bin"],wb);
      hdphiBJTall->Fill(m["dphiSL1"],wb);

      if (m["pairCodeSL1"]==0)
        hdphiBJTsig->Fill(m["dphiSL1"],wb);

      if (m["dphiSL1"]>PI23) {
        hmcxJAS->Fill(m["jtptSL"]/m["jtpt1"], m["pairCodeSL1"]==0 ? wb : wbkg);

        //signal xJ
        if (m["pairCodeSL1"]==0)
          hmcxJASsig->Fill(m["jtptSL"]/m["jtpt1"],wb);

        //signal-like pairCode (for MC purity)
        if (IsSignal(m))//m["subidSL"]==0)
          hPairCodeBFA->Fill(m["pairCodeSL1"],wb);////////////////////////////////there were w0 here
      }

      if (m["dphiSL1"]<PI13)
        hmcxJNS->Fill(m["jtptSL"]/m["jtpt1"],wb);


    }

//here I subtract background with TRUE leading jet!!!

float LJeff = 0.35;

if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSL"]>pt2cut
        && ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1)))
        hmcxJEars->Fill(m["jtptSL"]/m["jtpt1"],LJeff*(m["pairCodeSL1"]==0 ? wb : wbkg));

  });



  Fill(fdtinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2","hiHF"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["hiHF"]>5500) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCdata->Fill(m["dphi21"],w);

      if (dphi>PI23)
        hdtINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
      // if (dphi<PI13)
      //   hdtINCxJNS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
        hdtINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }


  });

  Fill(fmcinc,{"weight","jtpt1","refpt1","jtpt2","dphi21","bin","jteta1","jteta2","pairCodeSL1","discr_csvV1_1","jtptSL","dphiSL1","pthat","subidSL","refptSL","subid2"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
    if (m["pthat"]<65) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);




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




    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["jtpt2"]>pt2cut) {

      hdphiINCall->Fill(m["dphi21"],w);

      if (m["subid2"]==0)
        hdphiINCsig->Fill(m["dphi21"],w);

      if (dphi>PI23) {
          hmcINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);

          if (m["subid2"]==0)
            hmcINCxJASsig->Fill(m["jtpt2"]/m["jtpt1"],w);
      }
      
      if ((dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1))
          hmcINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50  && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>PI23) {
      if (m["pairCodeSL1"]<4 && IsSignal(m))//m["subidSL"]==0)
        hPairCodeQCD->Fill(m["pairCodeSL1"],w);
    }

  });


  hPairCode->SetBinContent(1,hPairCodeBFA->GetBinContent(1));
  hPairCode->SetBinContent(2,hPairCodeQCD->GetBinContent(2));
  hPairCode->SetBinContent(3,hPairCodeBFA->GetBinContent(3));
  hPairCode->SetBinContent(4,hPairCodeBFA->GetBinContent(4));
//  hPairCode->SetBinContent(5,hPairCodeBFA->GetBinContent(5)*30); //cheating, interesting idea

  // Normalize({hdtvz,hmcvz,hdtbin,hmcbin});

  // Draw({hdtvz,hmcvz});
  // plotylog = true;
  // Draw({hdtbin,hmcbin});

  //hdtxJASSub->Add(hdtxJAS,hdtxJNS,1,-1);
  //hmcxJASSub->Add(hmcxJAS,hmcxJNS,1,-1);

  //TODO: fix bin handling
  float coef = bkgfractionInNearSide[getbinindex((binMin+binMax)/2)];
  hdtxJASSubEars->Add(hdtxJAS,hdtxJEars,1,-1*coef);
  hmcxJASSubEars->Add(hmcxJAS,hmcxJEars,1,-1*coef);

  cout<<"   Fraction of data left after subtraction in bin "<<binMin/2<<" - "<<binMax/2<<" = "<<hdtxJASSubEars->Integral()/hdtxJAS->Integral()<<endl;

  hdtINCxJASSubEars->Add(hdtINCxJAS,hdtINCxJEars,1,-1);
  hmcINCxJASSubEars->Add(hmcINCxJAS,hmcINCxJEars,1,-1);

  //Ideal Hydjet+"Signal" subtraction:
  //vector<float> coef = {0.755140, 0.261050, 0.0134587};
  //calculated by simple dphi
  //vector<float> coef = {0.850222,0.435748,0.010241};
  //calculated by fancy
  //vector<float> coef = {0.869067,0.42076,0};

  //"quenched" away-side by 10%
  //vector<float> coef = {0.864643,0.490077,0.105538};

  for (int i=0;i<Nb;i++) {
    hsub[i]->Add(hasd[i],hbkg[i],1,-1*bkgfractionInNearSide[i]);
    hshj[i]->Add(hasd[i],hhyj[i],1,-1);
  }
  for (int i=0;i<Nbinc;i++) 
    hincsub[i]->Add(hincasd[i],hincbkg[i],1,-1);

  buildh(cbins);//Nb,0,100);
  auto hcentrSubSIG = geth("hcentrSubSIG","Signal;bin;<x_{J}>");
  auto hcentrSubASD = geth("hcentrSubASD","Unsubtracted;bin;<x_{J}>");
  auto hcentrSubBKG = geth("hcentrSubBKG","Background;bin;<x_{J}>");
  auto hcentrSubCLS = geth("hcentrSubCLS","Subtracted;bin;<x_{J}>");
  auto hcentrSubHJS = geth("hcentrSubHJS","Subtracted Hydjet;bin;<x_{J}>");


  buildh(Nbinc,0,100);
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
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  
  plotytitle = "Event fractions";

  plottextposx = 10;
  plotlegendpos = TopRight;
  plotputmean = false;
  plotymin = 0.;

  for (int i=0;i<Nb;i++) {
          plotymax = 1E-8;
    Draw({hbkg[i],hhyj[i]});
          plotymax = 1E-7;
    Draw({hsig[i],hasd[i],hsub[i],hshj[i]});
  }

for (int i=0;i<Nb;i++)
  cout<<"Amount of background in the subtraction "<<binnames[i]<<" : "<<
      setprecision(2)<<(1-hptNSsig[i]->Integral(0,hptNSsig[i]->GetNbinsX()+1)/(float)hptNS[i]->Integral(0,hptNS[i]->GetNbinsX()+1))*100<<" % "<<endl;;



plotymax = 0.5;//2.5E-8;
Normalize({hptNSsig[0],hptNSsig[1],hptNSsig[2],hptNS[0],hptNS[1],hptNS[2]});
  Draw({hptNSsig[0],hptNSsig[1],hptNSsig[2]});
  Draw({hptNS[0],hptNS[1],hptNS[2]});


  plotymin = 0.5;//0.4;
  plotymax = 0.8;//0.8;
  plotlegendpos = BottomRight;
  plottextposx = 0.5;
  plottextposy = 0.79;

  plotputmean = false;
  Draw({hcentrSubSIG, hcentrSubASD, hcentrSubCLS,hcentrSubHJS});



  Draw({hcentrSubSIGinc, hcentrSubASDinc, hcentrSubCLSinc});



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

  Draw({hmcINCxJAS,hmcINCxJEars,hmcINCxJASSubEars,hmcINCxJASsig});
  Draw({hmcxJAS,hmcxJEars,hmcxJASSubEars,hmcxJASsig});

  DrawCompare(hmcINCxJAS,hmcINCxJEars);
  DrawCompare(hmcINCxJASSubEars,hmcINCxJASsig);
  DrawCompare(hmcxJAS,hmcxJEars);
  DrawCompare(hmcxJASSubEars,hmcxJASsig);




//oops
  
  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcINCxJASSubEars});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsig});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsig});





  NormalizeAllHists();


  plotputmean = false;
  plotputwidth = true;
    plotymax = 0.4;

  DrawCompare(hdphiINCdata,hdphiINCall);
  DrawCompare(hdphiINCdata,hdphiINCsig);

  DrawCompare(hdphiBJTdata,hdphiBJTall);
  DrawCompare(hdphiBJTdata,hdphiBJTsig);

  Draw({hdphiINCdata,hdphiINCall,hdphiINCsig});
  Draw({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});
  plotputmean = true;
  plotputwidth = false;

  // plotputmean = true;
  // plotymax = 0.26;

  // Draw({hmcINCxJASsig,hmcINCxJEars,hmcINCxJAS,hmcINCxJASSubEars});
  // Draw({hmcxJASsig,hmcxJEars,hmcxJAS,hmcxJASSubEars});



  DrawCompare(hdtxJASSubEars, hmcxJASSubEars);

  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %%",binMin/2, binMax/2);

  DrawCompare(hdtINCxJASSubEars,hmcINCxJASSubEars);

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
  xjmcbjtmean.push_back(       hmcxJASSubEars->GetMean());
  xjmcbjtmeanerror.push_back(  hmcxJASSubEars->GetMeanError());

}



void tellmetruth()
{
  processWeights[0] = 1.2;
  processWeights[1] = 1.;
  processWeights[2] = 0.04;
  processWeights[3] = 0.04;

  //loadTagEffCorrections();


  findtruthpp();
  findtruthPbPb(60,200);
  findtruthPbPb(20,60); 
  findtruthPbPb(0,20);
  findtruthPbPb(0,200);

    //findtruthPbPb(140,200);
  
  cout<<"Bin \t\tInc.Data   \tInc.MC    \tSL.Data   \tSL.MC"<<endl;

  for (unsigned i=0;i<lbinmin.size();i++)
    cout<<setprecision(3)<<(int)lbinmin[i]/2<<" - "<<(int)lbinmax[i]/2<<" : \t"<<xjdtincmean[i]<<"\t"<<xjdtincmeanerror[i]<<
                                          "\t"<<xjmcincmean[i]<<"\t"<<xjmcincmeanerror[i]<<
                                          "\t"<<xjdtbjtmean[i]<<"\t"<<xjdtbjtmeanerror[i]<<
                                          "\t"<<xjmcbjtmean[i]<<"\t"<<xjmcbjtmeanerror[i]<<endl;
  

}
