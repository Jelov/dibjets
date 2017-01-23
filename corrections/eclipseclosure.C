#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "eclipsecorrections.h"

bool bjets = true;

void Plot(TString filename, TH1F *hxjAStrue,TH1F *hxjASsiguncorr, TH1F *hxjASsigcorr)
{

  float xjtrue = hxjAStrue->GetMean();
  float xjuncorr = hxjASsiguncorr->GetMean();
  float xjcorr = hxjASsigcorr->GetMean();

  float exjtrue = hxjAStrue->GetMeanError();
  float exjuncorr = hxjASsiguncorr->GetMeanError();
  float exjcorr = hxjASsigcorr->GetMeanError();

  auto c = getc();
  hxjAStrue->SetMinimum(0);
  hxjAStrue->SetMaximum(0.4);
  hxjAStrue->SetLineWidth(2);
  hxjAStrue->SetMarkerStyle(kNone);
  hxjAStrue->SetFillStyle(0);

  hxjASsiguncorr->SetMarkerStyle(kOpenSquare);
  hxjASsigcorr->SetMarkerStyle(kFullCircle);
  if (bjets)
    SetB({hxjAStrue,hxjASsiguncorr,hxjASsigcorr});
  else
    SetInc({hxjAStrue,hxjASsiguncorr,hxjASsigcorr});


  hxjAStrue->Draw("hist");
  hxjASsiguncorr->Draw("E1,same");
  hxjASsigcorr->Draw("E1,same");

  plotlegendpos = TopLeft;
  auto l = getLegend();
  l->AddEntry(hxjAStrue,Form("%s #LTx_{J}#GT=%.3f#pm%.3f",hxjAStrue->GetTitle(),xjtrue,exjtrue),"L");
  l->AddEntry(hxjASsiguncorr,Form("%s #LTx_{J}#GT=%.3f#pm%.3f",hxjASsiguncorr->GetTitle(),xjuncorr,exjuncorr),"P");
  l->AddEntry(hxjASsigcorr,Form("%s #LTx_{J}#GT=%.3f#pm%.3f",hxjASsigcorr->GetTitle(),xjcorr,exjcorr),"P");
  l->Draw();
  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.2, 0.8, aktstring);
  SavePlot(c,filename);

}

vector<float> xjtrue, xjsubuncorr,xjunsub,xjsubcorr,xjunsubcorr;
vector<float> exjtrue, exjsubuncorr,exjunsub,exjsubcorr,exjunsubcorr;

void checkclosure(int binMin, int binMax)
{

  auto f = config.getfile_djt(bjets ? "mcPbbfa" :"mcPbqcd");

  seth(10,0,1);
  auto hxjASsiguncorr = geth("hxjASsiguncorr","Signal, not eclipsed;x_{J};Event fractions");
  auto hxjASsigcorr = geth("hxjASsigcorr","Signal, corrected;x_{J};Event fractions");

  auto hxjAStrue = geth("hxjAStrue","Signal, full;x_{J};Event fractions");
  auto hxjAScorr = geth("hxjAScorr","Away-side weighted;x_{J};Event fractions");
  auto hxjNScorr = geth("hxjNScorr","Near-side weighted;x_{J};Event fractions");

  auto hxjASuncorr = geth("hxjASuncorr","Before subtraction;x_{J};Event fractions");//Away-side uncorrected
  auto hxjNSuncorr = geth("hxjNSuncorr","Near-side uncorrected;x_{J};Event fractions");
  auto hxjbkgsubuncorr = geth("hxjbkgsubuncorr","Bkg subtracted, uncorrected;x_{J};Event fractions");//

  auto hxjASclos = geth("hxjASclos","Bkg subtracted, corrected;x_{J};Event fractions");

  seth(1000,0,1000);
  auto hcorr = geth("hcorr");
  seth(1000,1,100,200,0,200);
  auto hcorr2d = geth2d("hcorr2d");

// nt->Draw("jtpt2/jtpt1>>h0(10,0,1)","weight*(jtpt1>100 && jtpt2>40 && bin<20 && pthat>65 && subid2==0 && refpt2>20)")
// nt->Draw("jtptSignal2/jtpt1>>h(10,0,1)","weight*(jtpt1>100 && jtptSignal2>40 && bin<20 && pthat>65 && dphiSignal21>2.1)")
// nt->Draw("jtpt2/jtpt1>>h3(10,0,1)","weight*(jtpt1>100 && jtpt2>40 && bin<20 && pthat>65 && dphi21>2.1)/(9.99514e-01*exp(-1.62510e+02*exp(-9.69483e-02*jtpt2)))")
// nt->Draw("jtpt2/jtpt1>>h2(10,0,1)","weight*(jtpt1>100 && jtpt2>40 && bin<20 && pthat>65 && dphi21<1.05)/(9.99514e-01*exp(-1.62510e+02*exp(-9.69483e-02*jtpt2)))")


  unordered_set<int> eventstomiss = {769779,1551232,1805770,1116573,1084397};

  Fill(f, [&](dict &d) {
    if (d["pthat"]<pthatcut) return;
    float w = d["weight"];
    float wsig = d["weight"];

    if (bjets && eventstomiss.find((int)d["event"])!=eventstomiss.end()) return;

    float corr = eclipseWeightmc(d["jtpt2"],d["bin"]);
    float nsfraction = NSfracmc(d["bin"]);//bjets ? 1 : 

    hcorr->Fill(eclipseWeightmc(d["jtpt2"],d["bin"]),w);
    hcorr2d->Fill(eclipseWeightmc(d["jtpt2"],d["bin"]),d["bin"],w);
    if (d["bin"]<binMin || d["bin"]>=binMax) return;
    bool dijetneeded = !bjets || (d["discr_csvV1_1"]>0.9 && d["discr_csvV1_2"]>0.9); // d["pairCode21"]==0;//
    bool signaldijetneeded = !bjets || (d["discr_csvV1_1"]>0.9 && d["discr_csvV1_Signal2"]>0.9); //d["pairCodeSignal21"]==0;//

    if (bjets) {
      if (d["pairCode21"]==0) w*=processweight(d["bProdCode"]);
      if (d["pairCodeSignal21"]==0) wsig*=processweight(d["bProdCode"]);
    }

    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && dijetneeded) {
      if (d["subid2"]==0 && d["refpt2"]>20 && d["dphi21"]>PI23) //
      {
        hxjASsiguncorr->Fill(d["jtpt2"]/d["jtpt1"],w);
        hxjASsigcorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr);
      }

      if (d["dphi21"]>PI23) {
        hxjASuncorr->Fill(d["jtpt2"]/d["jtpt1"],w);
        hxjAScorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr);
      }

      if (d["dphi21"]<PI13){ // && !(d["subid2"]==0 && d["refpt2"]>20) //for bkg only
        hxjNSuncorr->Fill(d["jtpt2"]/d["jtpt1"],w*nsfraction);
        hxjNScorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr*nsfraction);
      }
    }

    if (d["jtpt1"]>pt1cut && d["jtptSignal2"]>pt2cut &&  d["dphiSignal21"]>PI23 && signaldijetneeded)
      hxjAStrue->Fill(d["jtptSignal2"]/d["jtpt1"],wsig);
  });
 
hxjASclos->Add(hxjAScorr,hxjNScorr,1,-1); //subtraction is done by weights
hxjbkgsubuncorr->Add(hxjASuncorr,hxjNSuncorr,1,-1);

NormalizeAllHists();

plotylog = true;
Draw({hcorr});
plotylog = false;

auto c = getc();
hcorr2d->Draw("colz");
c->SetLogz();
c->SetLogx();
SavePlot(c,"hcorr2d");


plotlegendpos = TopLeft;
plotputmean = true;
plotymax = 0.35;
  aktstring = TString::Format("PbPb %d-%d %%",binMin/2, binMax/2);


Draw({hxjASsiguncorr,hxjAStrue,hxjAScorr,hxjNScorr,hxjASclos});

Draw({hxjAScorr,hxjNScorr,hxjASclos});

Draw({hxjASsiguncorr,hxjAScorr,hxjASclos});

Draw({hxjASsiguncorr,hxjAStrue,hxjASclos});



Plot(Form("closuresignal%d%d",binMin,binMax),hxjAStrue,hxjASsiguncorr,hxjASsigcorr);
Plot(Form("closurewithbkgsub%d%d",binMin,binMax),hxjAStrue,hxjbkgsubuncorr,hxjASclos);
hxjbkgsubuncorr->SetTitle("After subtraction");
Plot(Form("closureofbkgsubonly%d%d",binMin,binMax),hxjASsiguncorr,hxjASuncorr,hxjbkgsubuncorr);


xjtrue.push_back(hxjAStrue->GetMean());
exjtrue.push_back(hxjAStrue->GetMeanError());

xjsubuncorr.push_back(hxjbkgsubuncorr->GetMean());
exjsubuncorr.push_back(hxjbkgsubuncorr->GetMeanError());

xjunsub.push_back(hxjASuncorr->GetMean());
exjunsub.push_back(hxjASuncorr->GetMeanError());

xjsubcorr.push_back(hxjASclos->GetMean());
exjsubcorr.push_back(hxjASclos->GetMeanError());

xjunsubcorr.push_back(hxjAScorr->GetMean());
exjunsubcorr.push_back(hxjAScorr->GetMeanError());

}

void fillhistbins(TH1F *h, vector<float> x, vector<float> e)
{
  for (unsigned i=0;i<x.size();i++) {
    h->SetBinContent(i+1,x[i]);
    h->SetBinError(i+1,e[i]);
  }
}

void drawallcurves()
{
  loadeclipse();

  auto c = getc();
  plotlegendpos = BottomRight;
  auto l = new TLegend(0.66,0.3,0.83,0.7);
  gStyle->SetPalette(kRainBow);

  for (unsigned i=0;i<binbounds.size()-1;i++) {
    auto  fecl = new TF1(Form("feclipse%d",i),"exp(-[0]*exp(-[1]*x))",40,140);
    fecl->SetParameters(par0dt[i],par1dt[i]);

    // fecl->SetLineColorAlpha(kgreen,1-0.9*((float)i)/binbounds.size());
    fecl->SetLineColor(TColor::GetColorPalette(i*30));
    l->AddEntry(fecl, FloatToStr(binbounds[i]/2.)+"-"+FloatToStr(binbounds[i+1]/2.)+"%","L");
    fecl->SetMaximum(1);
    fecl->Draw(i==0 ? "" : "same");

    fecl->GetXaxis()->SetTitle("p_{T} [GeV]");
    fecl->GetYaxis()->SetTitle("Subleading jet findiding efficiency");
    fecl->GetXaxis()->CenterTitle();
    fecl->GetYaxis()->CenterTitle();
  }
  l->Draw();

  CMS_lumi(c, iPeriod, 33 ); 

  SavePlot(c,"eclipsecurves");

}

void eclipseclosure(int eclmode = 0, int bkgsubtractionmode = 0, bool incjetBJET = false)
{
  eclipsemode = eclmode;
  bkgsubmode = bkgsubtractionmode;
  bjets = incjetBJET;

  TString incbjt = bjets ? "bjt":"inc";

  macro m(Form("eclipseclosure0830_%d_%d_%s",eclipsemode,bkgsubmode,incbjt.Data()));

  drawallcurves();

  checkclosure(0,20);
  checkclosure(20,60);
  checkclosure(60,200);

  seth(3,0,3);
  auto hxjtrue = geth("hxjtrue","Truth");
  auto hxjunsub = geth("hxjunsub","Raw");
  auto hxjsubuncorr = geth("hxjsubuncorr","Sideband subtracted");
  auto hxjsubcorr = geth("hxjsubcorr","+ Eclipse correction");
  auto hxjunsubcorr = geth("hxjunsubcorr","Raw + Eclipse correction");

  fillhistbins(hxjtrue, xjtrue, exjtrue);
  fillhistbins(hxjsubuncorr, xjsubuncorr, exjsubuncorr);
  fillhistbins(hxjunsub, xjunsub, exjunsub);
  fillhistbins(hxjsubcorr, xjsubcorr, exjsubcorr);
  fillhistbins(hxjunsubcorr,xjunsubcorr,exjunsubcorr);


  vector<TString> axnames;
 for (auto s:binnames) axnames.push_back(s);
  RenameBinLabelsX(hxjtrue,axnames);
  RenameBinLabelsX(hxjsubuncorr,axnames);
  RenameBinLabelsX(hxjunsub,axnames);
  RenameBinLabelsX(hxjsubcorr,axnames);  
  RenameBinLabelsX(hxjunsubcorr,axnames);  



  plotlegendpos = TopLeft;
  // Draw(hmcPbdphi);

  aktstring = "";

  ShuffleBins(hxjtrue,{3,2,1});
  ShuffleBins(hxjsubuncorr,{3,2,1});
  ShuffleBins(hxjunsub,{3,2,1});
  ShuffleBins(hxjsubcorr,{3,2,1});
  ShuffleBins(hxjunsubcorr,{3,2,1});

  plotymin = 0.5;
  plotymax = 0.8;

  plotputmean = false;

if (bjets)
  SetB({hxjtrue,hxjsubuncorr,hxjunsub,hxjsubcorr,hxjunsubcorr});
else
  SetInc({hxjtrue,hxjsubuncorr,hxjunsub,hxjsubcorr,hxjunsubcorr});

hxjtrue->SetMarkerColor(kBlack); hxjtrue->SetLineColor(kBlack);
hxjsubuncorr->SetMarkerColor(kBlack); hxjsubuncorr->SetLineColor(kBlack);

hxjtrue->SetMarkerSize(2.);
hxjsubuncorr->SetMarkerSize(1.4);
hxjunsub->SetMarkerSize(1.4);
hxjsubcorr->SetMarkerSize(1.4);
hxjunsubcorr->SetMarkerSize(1.4);


hxjtrue->SetMarkerStyle(kOpenDiamond);
hxjsubuncorr->SetMarkerStyle(kOpenSquare);
hxjunsub->SetMarkerStyle(kOpenCircle);
hxjsubcorr->SetMarkerStyle(kFullCircle);
hxjunsubcorr->SetMarkerStyle(kOpenSquare);

  hxjsubcorr->SetMinimum(0.6);
  hxjsubcorr->SetMaximum(0.75);

  auto c=getc();
  auto l = getLegend();
  l->SetHeader(bjets ? "b-jets" : "Inclusive jets");
  l->AddEntry(hxjtrue,hxjtrue->GetTitle(),"P");
  l->AddEntry(hxjunsub,hxjunsub->GetTitle(),"P");
  l->AddEntry(hxjsubuncorr,hxjsubuncorr->GetTitle(),"P");
  l->AddEntry(hxjsubcorr,hxjsubcorr->GetTitle(),"P");
  // l->AddEntry(hxjunsubcorr,hxjunsubcorr->GetTitle(),"P");


  hxjsubcorr->Draw();
  hxjtrue->Draw("same");
  hxjsubuncorr->Draw("same");
// hxjunsubcorr->Draw("same");
  hxjunsub->Draw("same");
  l->Draw();


  SavePlot(c,"eclipsedemonstration");
  // Draw({hxjtrue,hxjsubuncorr,hxjunsub,hxjsubcorr});

}