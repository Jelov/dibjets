#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "tageffcorrections.h"


void findtruthpp()
{
  buildh(10,0,1);
  auto hmcppxJTrue = geth("hmcppxJTrue","true;x_{J};Event fractions");
  auto hmcppxJTrueTag = geth("hmcppxJTrueTag","true tagged;x_{J};Event fractions");
  auto hmcppxJTrueTagCorr = geth("hmcppxJTrueTagCorr","true tagged corrected;x_{J};Event fractions");


  TFile *fmcpp = new TFile(config.getFileName_djt("mcppbfa"));
  Fill(fmcpp,{"pthat","weight","jtpt1","jteta1","refpt1","discr_csvV1_1","jtptSB","jtetaSB","dphiSB1","bProdCode", "discr_csvV1_SB","refparton_flavorForB1","pairCodeSB1"},[&] (dict &m) {
    if (m["pthat"]<pthatcut) return;
    // if (m["pthat"]<80) return;
    
    // if (m["bProdCode"]!=1) return;
    float w = m["weight"]*processweight((int)m["bProdCode"]);

    // float w = m["weight"];
    // if (m["bProdCode"]==2) return;

    float corr = getppcorrection(m["jtpt1"],m["jteta1"],m["jtptSB"],m["jtetaSB"]);
    float wb = w*corr;

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23)
      hmcppxJTrue->Fill(m["jtptSB"]/m["jtpt1"],w);

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["dphiSB1"]>PI23
        && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_SB"]>0.9) {
      hmcppxJTrueTag->Fill(m["jtptSB"]/m["jtpt1"],w);
      hmcppxJTrueTagCorr->Fill(m["jtptSB"]/m["jtpt1"],wb);
    }


  });

  NormalizeAllHists();
  plotputmean = true;
  plotytitle = "Event fractions";
  plotdivide = false;
  // aktstring += "R=0.4 |#eta|<2.0";
  // plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  // plotthirdline = "#Delta#phi>2/3#pi";
  plottextposy = 0.8;
  plottextposx = 0.2;

  plotmeanposy = 0.43;
  plotymax = 0.2;
  plotymin = 0;
  plotlegendpos = BottomRight;//TopLeft;
  
  aktstring = "pp";
  SetMC({hmcppxJTrue,hmcppxJTrueTag,hmcppxJTrueTagCorr});
  SetData({hmcppxJTrue});
  Draw({hmcppxJTrue,hmcppxJTrueTag,hmcppxJTrueTagCorr});


  hmcppxJTrue->SetMinimum(0);
  hmcppxJTrue->SetMaximum(0.3);
  hmcppxJTrue->SetLineWidth(2);
  hmcppxJTrue->SetMarkerStyle(kNone);
  hmcppxJTrue->SetFillStyle(0);

  hmcppxJTrueTag->SetMarkerStyle(kOpenCircle);
  hmcppxJTrueTagCorr->SetMarkerStyle(kOpenSquare);


  plotymax = 0.3;

  SetB({hmcppxJTrue,hmcppxJTrueTag,hmcppxJTrueTagCorr});


  float xjtrue = hmcppxJTrue->GetMean();
  float xjtruetag = hmcppxJTrueTag->GetMean();
  float xjtruetagcorr = hmcppxJTrueTagCorr->GetMean();

  float exjtrue = hmcppxJTrue->GetMeanError();
  float exjtruetag = hmcppxJTrueTag->GetMeanError();
  float exjtruetagcorr = hmcppxJTrueTagCorr->GetMeanError();

  auto c = getc();
  hmcppxJTrue->Draw("hist");
  hmcppxJTrueTag->Draw("E1,same");
  hmcppxJTrueTagCorr->Draw("E1,same");

  plotlegendpos = TopLeft;
  auto l = getLegend();
  l->AddEntry(hmcppxJTrue,Form("b-dijets, #LTx_{J}#GT=%.3f#pm%.3f",xjtrue,exjtrue),"L");
  l->AddEntry(hmcppxJTrueTag,Form("uncorrected, #LTx_{J}#GT=%.3f#pm%.3f",xjtruetag,exjtruetag),"P");
  l->AddEntry(hmcppxJTrueTagCorr,Form("corrected, #LTx_{J}#GT=%.3f#pm%.3f",xjtruetagcorr,exjtruetagcorr),"P");
  l->Draw();
  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.2, 0.8, aktstring);
  SavePlots(c,"closurepp");
}


void findtruthPbPb(int binMin, int binMax)
{
  TFile *fmc = new TFile(config.getFileName_djt("mcPbbfa"));

  buildNamesuffix = TString::Format("_bin_%d_%d",binMin, binMax);
  //  buildTitlesuffix = TString::Format("%d-%d %%",binMin/2, binMax/2);

  buildh(10,0,1);
  auto hmcPbPbxJTrue = geth("hmcPbPbxJTrue","PbPb true;x_{J};Event fractions");
  auto hmcPbPbxJTrueTag = geth("hmcPbPbxJTrueTag","PbPb true tagged;x_{J};Event fractions");
  auto hmcPbPbxJTrueTagCorr = geth("hmcPbPbxJTrueTagCorr","PbPb true tagged corrected;x_{J};Event fractions");
  auto hmcPbPbxJTrueTagCorrPt = geth("hmcPbPbxJTrueTagCorrPt","PbPb true tagged corrected pt;x_{J};Event fractions");
  auto hmcPbPbxJTrueTagCorrEta = geth("hmcPbPbxJTrueTagCorrEta","PbPb true tagged corrected eta;x_{J};Event fractions");
  auto hmcPbPbxJTrueTagCorrBin = geth("hmcPbPbxJTrueTagCorrBin","PbPb true tagged corrected bin;x_{J};Event fractions");

  buildh(12,20,140);//10,40,100);
  auto hpt2true = geth("hpt2true","true;p_{T,2} GeV");
  auto hpt2truetag = geth("hpt2truetag","true tagged;p_{T,2} GeV");
  auto hpt2truetagovertrue = geth("hpt2truetagovertrue","true tagged/true;p_{T,2} GeV");
  auto hpt2truetagcorr = geth("hpt2truetagcorr","true tagged corrected;p_{T,2} GeV");
  auto hpt2truetagcorrovertrue = geth("hpt2truetagcorrovertrue","true tagged corrected/true;p_{T,2} GeV");

  buildh(10,100,200);
  auto hpt1true = geth("hpt1true","true;p_{T,1} GeV");
  auto hpt1truetag = geth("hpt1truetag","true tagged;p_{T,1} GeV");
  auto hpt1truetagovertrue = geth("hpt1truetagovertrue","true tagged/true;p_{T,1} GeV");
  auto hpt1truetagcorr = geth("hpt1truetagcorr","true tagged corrected;p_{T,1} GeV");
  auto hpt1truetagcorrovertrue = geth("hpt1truetagcorrovertrue","true tagged corrected/true;p_{T,1} GeV");

  buildh(10,0,200);
  auto hbintrue = geth("hbintrue","true;bin");
  auto hbintruetag = geth("hbintruetag","true tagged;bin");
  auto hbintruetagovertrue = geth("hbintruetagovertrue","true tagged/true;bin");
  auto hbintruetagcorr = geth("hbintruetagcorr","true tagged corrected;bin");
  auto hbintruetagcorrovertrue = geth("hbintruetagcorrovertrue","true tagged corrected/true;bin");

  buildh(10,-2,2);
  auto heta2true = geth("heta2true","true;#eta_{2}");
  auto heta2truetag = geth("heta2truetag","true tagged;#eta_{2}");
  auto heta2truetagovertrue = geth("heta2truetagovertrue","true tagged/true;#eta_{2}");
  auto heta2truetagcorr = geth("heta2truetagcorr","true tagged corrected;#eta_{2}");
  auto heta2truetagcorrovertrue = geth("heta2truetagcorrovertrue","true tagged corrected/true;#eta_{2}");

  auto heta1true = geth("heta1true","true;#eta_{1}");
  auto heta1truetag = geth("heta1truetag","true tagged;#eta_{1}");
  auto heta1truetagovertrue = geth("heta1truetagovertrue","true tagged/true;#eta_{1}");
  auto heta1truetagcorr = geth("heta1truetagcorr","true tagged corrected;#eta_{1}");
  auto heta1truetagcorrovertrue = geth("heta1truetagcorrovertrue","true tagged corrected/true;#eta_{1}");

  Fill(fmc,{"bin","pthat","weight","jtpt1","jteta1","refpt1","discr_csvV1_1","jtptSB","refptSB","jtetaSB","dphiSB1","bProdCode", "discr_csvV1_SB","refparton_flavorForB1","pairCodeSB1"},[&] (dict &m) {
    if (m["bin"]<binMin || m["bin"]>=binMax) return;
    if (m["pthat"]<pthatcut) return;
    // if (m["pthat"]<80) return;

    // if (m["bProdCode"]!=1) return;
    float w = m["weight"]*processweight((int)m["bProdCode"]); //because we have only leading B and second B

    // float w = m["weight"];
    // if (m["bProdCode"]==2) return;


    float corr = getPbPbcorrection(m["jtpt1"],m["jteta1"],m["jtptSB"],m["jtetaSB"],m["bin"]);
    // float corrpt  = getPbPbcorrectionPt(m["jtpt1"],m["jtptSB"]);
    // float correta = getPbPbcorrectionEta(m["jteta1"],m["jtetaSB"]);
    // float corrbin = getPbPbcorrectionBin(m["bin"]);


    float wb = w*corr;



   if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["refptSB"]>20 && m["dphiSB1"]>PI23) {
      hmcPbPbxJTrue->Fill(m["jtptSB"]/m["jtpt1"],w);
      hpt2true->Fill(m["jtptSB"],w);
      hpt1true->Fill(m["jtpt1"],w);
      heta2true->Fill(m["jtetaSB"],w);
      heta1true->Fill(m["jteta1"],w);
      hbintrue->Fill(m["bin"],w);
   }

    if (m["jtpt1"]>pt1cut && m["refpt1"]>50 && abs(m["refparton_flavorForB1"])==5 && m["jtptSB"]>pt2cut && m["refptSB"]>20 && m["dphiSB1"]>PI23
        && m["discr_csvV1_1"]>0.9 &&  m["discr_csvV1_SB"]>0.9) { //

     //corrpt *= m["jtptSB"] < 60 ? 1./0.7 : 1;
     //wb *= m["jtptSB"] < 60 ? 1./0.7 : 1;

      hmcPbPbxJTrueTag->Fill(m["jtptSB"]/m["jtpt1"],w);
      hmcPbPbxJTrueTagCorr->Fill(m["jtptSB"]/m["jtpt1"],wb);
      // hmcPbPbxJTrueTagCorrPt->Fill(m["jtptSB"]/m["jtpt1"],w*corrpt);
      // hmcPbPbxJTrueTagCorrEta->Fill(m["jtptSB"]/m["jtpt1"],w*corrpt*correta);
      // hmcPbPbxJTrueTagCorrBin->Fill(m["jtptSB"]/m["jtpt1"],w*corrpt*correta*corrbin);


      hpt2truetag->Fill(m["jtptSB"],w);
      hpt1truetag->Fill(m["jtpt1"],w);
      heta2truetag->Fill(m["jtetaSB"],w);
      heta1truetag->Fill(m["jteta1"],w);
      hbintruetag->Fill(m["bin"],w);

      hpt2truetagcorr->Fill(m["jtptSB"],wb);
      hpt1truetagcorr->Fill(m["jtpt1"],wb);
      heta2truetagcorr->Fill(m["jtetaSB"],wb);
      heta1truetagcorr->Fill(m["jteta1"],wb);
      hbintruetagcorr->Fill(m["bin"],wb);

    }






  });

  NormalizeAllHists();
//plotymax = 9999;
  aktstring = TString::Format("PbPb %d-%d %%",binMin/2, binMax/2);//TString::Format("PbPb#Delta#phi>2/3#pi %d-%d %%",binMin/2, binMax/2);

  SetMC({hmcPbPbxJTrue,hmcPbPbxJTrueTag,hmcPbPbxJTrueTagCorr});
  SetData({hmcPbPbxJTrue});
  hmcPbPbxJTrue->SetMinimum(0);
  hmcPbPbxJTrue->SetMaximum(0.3);
  hmcPbPbxJTrue->SetLineWidth(2);
  hmcPbPbxJTrue->SetMarkerStyle(kNone);
  hmcPbPbxJTrue->SetFillStyle(0);

  hmcPbPbxJTrueTag->SetMarkerStyle(kOpenCircle);
  hmcPbPbxJTrueTagCorr->SetMarkerStyle(kOpenSquare);


  plotymax = 0.3;
  Draw({hmcPbPbxJTrue,hmcPbPbxJTrueTag,hmcPbPbxJTrueTagCorr});

  SetB({hmcPbPbxJTrue,hmcPbPbxJTrueTag,hmcPbPbxJTrueTagCorr});

  float xjtrue = hmcPbPbxJTrue->GetMean();
  float xjtruetag = hmcPbPbxJTrueTag->GetMean();
  float xjtruetagcorr = hmcPbPbxJTrueTagCorr->GetMean();

  float exjtrue = hmcPbPbxJTrue->GetMeanError();
  float exjtruetag = hmcPbPbxJTrueTag->GetMeanError();
  float exjtruetagcorr = hmcPbPbxJTrueTagCorr->GetMeanError();

  auto c = getc();
  hmcPbPbxJTrue->Draw("hist");
  hmcPbPbxJTrueTag->Draw("E1,same");
  hmcPbPbxJTrueTagCorr->Draw("E1,same");

  plotlegendpos = TopLeft;
  auto l = getLegend();
  l->AddEntry(hmcPbPbxJTrue,Form("b-dijets, #LTx_{J}#GT=%.3f#pm%.3f",xjtrue,exjtrue),"L");
  l->AddEntry(hmcPbPbxJTrueTag,Form("uncorrected, #LTx_{J}#GT=%.3f#pm%.3f",xjtruetag,exjtruetag),"P");
  l->AddEntry(hmcPbPbxJTrueTagCorr,Form("corrected, #LTx_{J}#GT=%.3f#pm%.3f",xjtruetagcorr,exjtruetagcorr),"P");
  l->Draw();
  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.2, 0.8, aktstring);
  SavePlots(c,Form("closure%d%d",binMin,binMax));


    // //if (binMin==0 && binMax==200) {

    // Draw({hmcPbPbxJTrueTag,hmcPbPbxJTrueTagCorrPt,hmcPbPbxJTrueTagCorrEta,hmcPbPbxJTrueTagCorrBin});
 

    // SetMC({hpt2truetag,hpt1truetag,heta2truetag,heta1truetag,hbintruetag});

    // plotputmean = false;

    // plotymax = 0.2;

    // Draw({hpt2true,hpt2truetag,hpt2truetagcorr});

    // plotymax = 0.3;

    // Draw({hpt1true,hpt1truetag,hpt1truetagcorr});

    // plotymax = 0.2;
    // Draw({heta2true,heta2truetag,heta2truetagcorr});
    // Draw({heta1true,heta1truetag,heta1truetagcorr});

    // plotymax = 1;
    // Draw({hbintrue,hbintruetag,hbintruetagcorr});


plotymin = 0;
plotymax = 0.2;

hpt2truetagovertrue->Divide(hpt2truetag,hpt2true,1,1,"B");
hpt1truetagovertrue->Divide(hpt1truetag,hpt1true,1,1,"B");
heta2truetagovertrue->Divide(heta2truetag,heta2true,1,1,"B");
heta1truetagovertrue->Divide(heta1truetag,heta1true,1,1,"B");
hbintruetagovertrue->Divide(hbintruetag,hbintrue,1,1,"B");


hpt2truetagcorrovertrue->Divide(hpt2truetagcorr,hpt2true,1,1,"B");
hpt1truetagcorrovertrue->Divide(hpt1truetagcorr,hpt1true,1,1,"B");
heta2truetagcorrovertrue->Divide(heta2truetagcorr,heta2true,1,1,"B");
heta1truetagcorrovertrue->Divide(heta1truetagcorr,heta1true,1,1,"B");
hbintruetagcorrovertrue->Divide(hbintruetagcorr,hbintrue,1,1,"B");

 NormalizeAllHists();

Draw({hpt2truetagovertrue,hpt2truetagcorrovertrue});
Draw({hpt1truetagovertrue,hpt1truetagcorrovertrue});
Draw({heta2truetagovertrue,heta2truetagcorrovertrue});
Draw({heta1truetagovertrue,heta1truetagcorrovertrue});
Draw({hbintruetagovertrue,hbintruetagcorrovertrue});



 // }

}



void tageffcorr()
{

  macro m("taggeffcorr_20160512");
  //looptupledryrun = false;
  //applyCorrection = true;


  loadTagEffCorrections();

  findtruthpp();
  findtruthPbPb(60,200);
  findtruthPbPb(20,60); 
  findtruthPbPb(0,20);
  findtruthPbPb(0,200);


}
