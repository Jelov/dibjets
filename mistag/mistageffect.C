#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"

const float pt1cut = 100;
const float pt2cut = 40;

const float pi23 = 3.142*2/3;
const float pi13 = 3.142*1/3;

vector<float> processWeights(4);

void effect(TString prefix, TString bfafile, TString qcdfile)
{
  buildh(10,0,1);
  auto hmcppxJTrue = geth(prefix+"hmcppxJTrue","true;x_{J}");
  auto hmcppxJMistag = geth(prefix+"hmcppxJMistag","mistag;x_{J}");
  auto hmcppxJTotal = geth(prefix+"hmcppxJTotal","total;x_{J}");
  buildh(5,0,5);
  auto puritybfa = geth(prefix+"puritybfa");
  auto purityqcd = geth(prefix+"purityqcd");
  auto purity = geth(prefix+"purity");



  TFile *fmcpp = new TFile(config.getFileName_djt(bfafile));
  Fill(fmcpp,{"bProdCode","pthat","weight","jtpt1","refpt1","discr_csvV1_1","refparton_flavorForB1","jtptSL","subidSL","refparton_flavorForBSL","pairCodeSL1","dphiSL1",},[&] (dict &m) {
    if (m["pthat"]<50) return;
    if ((abs(m["refparton_flavorForB1"])!=5 && abs(m["refparton_flavorForBSL"])!=5)|| m["subidSL"]!=0 ) return;// || m["subidSL"]!=0

    float w = m["weight"];
    float wb = w*processWeights[(int)m["bProdCode"]];

    //    w = wb;
    //float corr = getppcorrection(m["jtpt1"],m["jteta1"],m["jtptSB"],m["jtetaSB"]);
    //float wb = w*corr;

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23 && m["discr_csvV1_1"]>0.9) { //discr of SL is always >0.9
      hmcppxJTotal->Fill(m["jtptSL"]/m["jtpt1"],wb);

      if (m["pairCodeSL1"]==0) {
        hmcppxJTrue->Fill(m["jtptSL"]/m["jtpt1"],wb);
	puritybfa->Fill(m["pairCodeSL1"],wb);
      }
      else{
        hmcppxJMistag->Fill(m["jtptSL"]/m["jtpt1"],w);
	puritybfa->Fill(m["pairCodeSL1"],w);
      }

    }

  });


  TFile *fmcppqcd = new TFile(config.getFileName_djt(qcdfile));
  Fill(fmcppqcd,{"bProdCode","pthat","weight","jtpt1","refpt1","discr_csvV1_1","refparton_flavorForB1","jtptSL","subidSL","refparton_flavorForBSL","pairCodeSL1","dphiSL1",},[&] (dict &m) {
    if (m["pthat"]<50) return;
    if (abs(m["refparton_flavorForB1"])==5 || abs(m["refparton_flavorForBSL"])==5 || m["pairCodeSL1"]==4 || m["subidSL"]!=0 ) return;// || m["subidSL"]!=0

    float w = m["weight"];
    //float corr = getppcorrection(m["jtpt1"],m["jteta1"],m["jtptSB"],m["jtetaSB"]);
    //float wb = w*corr;

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23 && m["discr_csvV1_1"]>0.9) { //discr of SL is always >0.9
      hmcppxJTotal->Fill(m["jtptSL"]/m["jtpt1"],w);

      hmcppxJMistag->Fill(m["jtptSL"]/m["jtpt1"],w);

      purityqcd->Fill(m["pairCodeSL1"],w);
    }

  });

  purity->Add(puritybfa,purityqcd);

  NormalizeAllHists();

  Print(purity);

  plotputmean = true;
  plotytitle = "Event fractions";
  plotdivide = false;
  // aktstring += "R=0.4 |#eta|<2.0";
  // plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  // plotthirdline = "#Delta#phi>2/3#pi";
  plottextposy = 0.8;
  plottextposx = 0.2;

  plotmeanposy = 0.43;
  plotymax = 0.3;
  plotymin = 0;
  plotlegendpos = BottomRight;//TopLeft;
  
  aktstring = prefix;
  SetMC({hmcppxJTrue,hmcppxJMistag,hmcppxJTotal});
  SetData({hmcppxJTrue});
  Draw({hmcppxJTrue,hmcppxJMistag,hmcppxJTotal});

  plotputmean = false;
  plotymax = 1.;
  Draw({purity},"bar");


}


void mistageffect()
{
  //looptupledryrun = false;
  //applyCorrection = true;

  processWeights[0] = 1.2;
  processWeights[1] = 1.;
  processWeights[2] = 0.04;
  processWeights[3] = 0.04;

  // processWeights[0] = 1.2;
  // processWeights[1] = 1.;
  // processWeights[2] = 0.04;
  // processWeights[3] = 0.04;

  // loadTagEffCorrections();

  effect("wFEXpp","mcppbfa","mcppqcd");
  effect("wFEXPbPb","mcPbbfa","mcPbqcd");

  // findtruthPbPb(60,200);
  // findtruthPbPb(20,60); 
  // findtruthPbPb(0,20);
  // findtruthPbPb(0,200);


}
