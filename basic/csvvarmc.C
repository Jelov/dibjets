#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "../corrections/tageffcorrections.h"
#include "../corrections/eclipsecorrections.h"

class listofhist
{
public:
  vector<TH1F *> hists;
  listofhist(vector<float> range, TString nameprefix)
  {

  }

};


vector<tageffcorr *> cor;
vector<float> csvrange = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1};
int ncsv = 0;
bool pp;
TString binstr;

void initcor(bool fix095)//(bool only0p9)
{


  // for (auto csv:csvrange)
  //   cor.push_back(new tageffcorr(only0p9 ? 0.9 : csv));

  cor.clear();
  for (unsigned i=0;i<csvrange.size();i++) {
    float c = csvrange[i];
    if (fix095 && i==csvrange.size()-2) c = 0.9;
    cout<<c<<endl;
    cor.push_back(new tageffcorr(c));
  }


}


void CumDefPlot(vector<float> numbers, TH1F *hist, TString name)
{
  vector<TGraph *> gg;
  auto c = getc();
  c->SetGrid();
  plotlegendpos = TopLeft;

  auto l = new TLegend(0.18,0.65,0.4,0.8);
  auto g = new TGraph(ncsv,&csvrange[0],&(numbers[0]));
  g->SetMinimum(0);
  g->SetMaximum(1);
  g->SetMarkerColor(kRed);
  gg.push_back(g);

  g->Draw("APL");
  g->GetXaxis()->SetTitle("CSV threshold");
  g->GetYaxis()->SetTitle(name);

  hist->Draw("same");

  l->SetHeader(binstr);
  l->AddEntry(g,"cumulative","P");
  l->AddEntry(hist,"differential","P");
  l->Draw();

  SavePlot(c,name+binstr);//(pp ? "_pp" : Form("_bin_%d_%d",binMin,binMax)));

}

TGraph *getgraph(TH1F *h)
{
  vector<float> x,y,dx,dy;
  for (int i=0;i<h->GetNbinsX();i++) {
    x.push_back(h->GetBinLowEdge(i+1));
    y.push_back(h->GetBinContent(i+1));
    dy.push_back(h->GetBinError(i+1));
    dx.push_back(0);
  }
  auto g = new TGraphErrors(x.size(),&x[0],&y[0],&dx[0],&dy[0]);
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetMarkerSize(h->GetMarkerSize());
  g->SetLineColor(h->GetLineColor());
  g->SetMinimum(h->GetMinimum());
  g->SetMaximum(h->GetMaximum());

  return g;

}

void makeplots(int binMin, int binMax, bool MC, float sys)
{
  pp = binMin==-1 && binMax==-1; //pp = false, PbPb = true
  
  binstr = pp ? "pp" : Form("%d-%d %%",binMin/2,binMax/2);

  TString filename;
  if (MC) filename = pp ? "mcppbfa":"mcPbbfa";
    else  filename = pp ? "dtppjpf":"dtPbbjt";

  auto fmcbfa = config.getfile_djt(filename);//pp ? "mcppbfa":"mcPbbfa");

  buildNamesuffix = pp ? "pp" : TString::Format("_bin_%d_%d",binMin, binMax);


  ncsv = csvrange.size();


  vector<TH1F *> vhPairCode,hxjcum,hxjcumcor,hxjdif;
  vector<TH1F *> hxjbkgcum,hxjbkgcumcor;
  vector<TH1F *> hxjallcum,hxjallcumcor;
  vector<TH1F *> hxjeclsigallcum,hxjeclsigallcumcor;
  vector<TH1F *> hxjeclallcum,hxjeclallcumcor;

  vector<TH1F *> hxjuneclsigallcum,hxjuneclsigallcumcor;
  vector<TH1F *> hxjuneclallcum,hxjuneclallcumcor;
  vector<TH1F *> hxjuneclallNScum,hxjuneclallNScumcor,hxjuneclallSUBcumcor;

  vector<float> effnum(ncsv);
  vector<float> combpurAS(ncsv),combpurNS(ncsv), combpurdenAS(ncsv),combpurdenNS(ncsv);
  float effden = 0;


  for (int i=0;i<ncsv;i++) {
    seth(5,0,5);

    vhPairCode.push_back(geth(Form("vhPairCode_%.2f",csvrange[i])));
  
    seth(10,0,1);/////////////////////////////////////////!!!!!!!!!!!!!!!!!!
    hxjcum.push_back(geth(Form("hxjcum%.2f",csvrange[i])));
    hxjcumcor.push_back(geth(Form("hxjcumcor%.2f",csvrange[i])));
    hxjdif.push_back(geth(Form("hxjdif%.2f",csvrange[i])));

    hxjbkgcum.push_back(geth(Form("hxjbkgcum%.2f",csvrange[i])));
    hxjbkgcumcor.push_back(geth(Form("hxjbkgcumcor%.2f",csvrange[i])));

    hxjallcum.push_back(geth(Form("hxjallcum%.2f",csvrange[i])));
    hxjallcumcor.push_back(geth(Form("hxjallcumcor%.2f",csvrange[i])));

hxjeclsigallcum.push_back(geth(Form("hxjeclsigallcum%.2f",csvrange[i])));
hxjeclsigallcumcor.push_back(geth(Form("hxjeclsigallcumcor%.2f",csvrange[i])));
hxjeclallcum.push_back(geth(Form("hxjeclallcum%.2f",csvrange[i])));
hxjeclallcumcor.push_back(geth(Form("hxjeclallcumcor%.2f",csvrange[i])));

hxjuneclsigallcum.push_back(geth(Form("hxjuneclsigallcum%.2f",csvrange[i])));
hxjuneclsigallcumcor.push_back(geth(Form("hxjuneclsigallcumcor%.2f",csvrange[i])));
hxjuneclallcum.push_back(geth(Form("hxjuneclallcum%.2f",csvrange[i])));
hxjuneclallcumcor.push_back(geth(Form("hxjuneclallcumcor%.2f",csvrange[i])));

hxjuneclallNScum.push_back(geth(Form("hxjuneclallNScum%.2f",csvrange[i])));
hxjuneclallNScumcor.push_back(geth(Form("hxjuneclallNScumcor%.2f",csvrange[i])));
hxjuneclallSUBcumcor.push_back(geth(Form("hxjuneclallSUBcumcor%.2f",csvrange[i])));


  combpurdenAS.push_back(0);
  combpurdenNS.push_back(0);
  }


  seth(csvrange);
  auto hpuritynum = geth("hpuritynum");
  auto hpurityden = geth("hpurityden");
  auto hpurity = geth("hpurity",";Subleading jet CSV threshold;purity");

  auto hcombpurityASnum = geth("hcombpurityASnum");
  auto hcombpurityASden = geth("hcombpurityASden");
  auto hcombpurityAS = geth("hcombpurityAS",";Subleading jet CSV threshold;combinatorial purity AS");

  auto hcombpurityNSnum = geth("hcombpurityNSnum");
  auto hcombpurityNSden = geth("hcombpurityNSden");
  auto hcombpurityNS = geth("hcombpurityNS",";Subleading jet CSV threshold;combinatorial purity NS");


  auto hefficiency = geth("hefficiency",";Subleading jet CSV threshold;efficiency");


  auto hmeanxjcum = geth("hmeanxjcum","cumulative;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjcumcor = geth("hmeanxjcumcor","cumulative corrected;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjdif = geth("hmeanxjdif","differential;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hmeanxjbkgcum = geth("hmeanxjbkgcum","bkg cumulative;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjbkgcumcor = geth("hmeanxjbkgcumcor","bkg cumulative corrected;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hmeanxjallcum = geth("hmeanxjallcum","all cumulative;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjallcumcor = geth("hmeanxjallcumcor","all cumulative corrected;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hxjtrue = geth("hxjtrue","true xJ");
  auto hxjeclipsedtrue = geth("hxjeclipsedtrue","eclipsed true xJ");


  auto hmeanxjeclsigallcum = geth("hmeanxjeclsigallcum","eclipsed signal all;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjeclsigallcumcor = geth("hmeanxjeclsigallcumcor","eclipsed signal all corrected;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hmeanxjeclallcum = geth("hmeanxjeclallcum","eclipsed signal all;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjeclallcumcor = geth("hmeanxjeclallcumcor","eclipsed signal all corrected;Subleading jet CSV threshold;#LTx_{J}#GT");


  auto hmeanxjuneclsigallcum = geth("hmeanxjuneclsigallcum","UN eclipsed signal all;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjuneclsigallcumcor = geth("hmeanxjuneclsigallcumcor","UN eclipsed signal all corrected;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hmeanxjuneclallcum = geth("hmeanxjuneclallcum","UN eclipsed signal all;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjuneclallcumcor = geth("hmeanxjuneclallcumcor","UN eclipsed signal all corrected;Subleading jet CSV threshold;#LTx_{J}#GT");

  auto hmeanxjuneclallNScum = geth("hmeanxjuneclallNScum","UN eclipsed all NOT corrected NS;Subleading jet CSV threshold;#LTx_{J}#GT");
  auto hmeanxjuneclallNScumcor = geth("hmeanxjuneclallNScumcor","UN eclipsed all corrected NS;Subleading jet CSV threshold;x_{J}");
  auto hmeanxjuneclallSUBcumcor = geth("hmeanxjuneclallSUBcumcor","UN eclipsed all corrected SUB;Subleading jet CSV threshold;#LTx_{J}#GT");

  Fill(fmcbfa,[&] (dict &d) {

    if (d["event"]==128751 || d["event"]==1551232 || d["event"]==3043866) return;

    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;
    if (MC && d["pthat"]<pthatcut) return;

    if (MC && abs(d["refparton_flavorForB1"])!=5) return;
    if (!MC && d["discr_csvV1_1"]<0.9) return;

    if (MC && d["refpt1"]<50) return;
    if (!(d["jtpt1"]>pt1cut)) return;

    float w = MC ? weight12(d) : d["weight"];
    if (d["numTagged"]>6) return;



//true signal2
    if (MC && d["jtptSignal2"]>pt2cut && d["dphiSignal21"]>PI23) {
      float ws2 = weight1Signal2(d);
      float csv2 = d["discr_csvV1_Signal2"];

      float xj = d["jtptSignal2"]/d["jtpt1"];
      float paircode = d["pairCodeSignal21"];

      if (paircode==0) {
        hxjtrue->Fill(xj,ws2);
        hefficiency->Fill(csv2,ws2);   
        effden+=ws2;  
        hpuritynum->Fill(csv2,ws2);   
      }
      hpurityden->Fill(csv2,ws2);
    
      for (int i=0;i<ncsv;i++)
        if (csv2>csvrange[i]) {

          float c = pp ? cor[i]->pp(d["jtpt1"],d["jteta1"],d["jtptSignal2"],d["jtetaSignal2"]):
                         cor[i]->PbPb(d["jtpt1"],d["jteta1"],d["jtptSignal2"],d["jtetaSignal2"],d["bin"]);

          vhPairCode[i]->Fill(paircode,ws2);

          if (paircode==0) {
            hxjcum[i]->Fill(xj,ws2);
            hxjcumcor[i]->Fill(xj,ws2*c);

            if (i==ncsv-1 || csv2<csvrange[i+1])
              hxjdif[i]->Fill(xj,ws2);

            effnum[i]+=ws2;
          }

          if (paircode!=0) {
            hxjbkgcum[i]->Fill(xj,ws2);
            hxjbkgcumcor[i]->Fill(xj,ws2*c);
          }

          hxjallcum[i]->Fill(xj,ws2);
          hxjallcumcor[i]->Fill(xj,ws2*c);

        }
    }

//away-side jets jt2

    float w2 = MC ? weight12(d) : d["weight"];
    float csv2 = d["discr_csvV1_2"];
    float ew = 1;
    if (!pp) ew = MC ? eclipseWeightmc(d["jtpt2"],d["bin"]) : eclipseWeightdt(d["jtpt2"],d["bin"]);
    float xj = d["jtpt2"]/d["jtpt1"];
    bool signaljet = MC && d["subid2"]==0 && d["refpt2"]>20;
    float paircode = MC ? d["pairCode21"] : -1;


    if (d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {



      if (signaljet && paircode==0)
        hxjeclipsedtrue->Fill(xj,w2);

      if (signaljet)
        hcombpurityASnum->Fill(csv2,w2);
      hcombpurityASden->Fill(csv2,w2);

      for (int i=0;i<ncsv;i++)
        if (csv2>csvrange[i]) {

          float c = pp ? cor[i]->pp(d["jtpt1"],d["jteta1"],d["jtpt2"],d["jteta2"]):
                         cor[i]->PbPb(d["jtpt1"],d["jteta1"],d["jtpt2"],d["jteta2"],d["bin"]);

          if (signaljet) {
            hxjeclsigallcum[i]->Fill(xj,w2);
            hxjeclsigallcumcor[i]->Fill(xj,w2*c);
            hxjuneclsigallcum[i]->Fill(xj,w2*ew);
            hxjuneclsigallcumcor[i]->Fill(xj,w2*c*ew);
            combpurAS[i]+=w2;

          }
          hxjeclallcum[i]->Fill(xj,w2);
          hxjeclallcumcor[i]->Fill(xj,w2*c);
          hxjuneclallcum[i]->Fill(xj,w2*ew);
          hxjuneclallcumcor[i]->Fill(xj,w2*c*ew);
          combpurdenAS[i]+=w2;

        }
    }

    if (d["jtpt2"]>pt2cut && d["dphi21"]<PI13) {

      if (signaljet)
        hcombpurityNSnum->Fill(csv2,w2);
      hcombpurityNSden->Fill(csv2,w2);

      for (int i=0;i<ncsv;i++)
        if (csv2>csvrange[i]) {
          float c = pp ? cor[i]->pp(d["jtpt1"],d["jteta1"],d["jtpt2"],d["jteta2"]):
                          cor[i]->PbPb(d["jtpt1"],d["jteta1"],d["jtpt2"],d["jteta2"],d["bin"]);
          hxjuneclallNScum[i]->Fill(xj,w2*ew);
          hxjuneclallNScumcor[i]->Fill(xj,w2*c*ew);

          if (signaljet) {
            combpurNS[i]+=w2;

          }
          combpurdenNS[i]+=w2;

        }
    }


  });

  hpurity->Divide(hpuritynum,hpurityden,1,1,"B");
  hcombpurityAS->Divide(hcombpurityASnum,hcombpurityASden,1,1,"B");  
  hcombpurityNS->Divide(hcombpurityNSnum,hcombpurityNSden,1,1,"B");  
  for (unsigned i=0;i<csvrange.size();i++)
  {
    hmeanxjcum->SetBinContent(i+1,hxjcum[i]->GetMean());
    hmeanxjcum->SetBinError(i+1,hxjcum[i]->GetMeanError());
    hmeanxjcumcor->SetBinContent(i+1,hxjcumcor[i]->GetMean());
    hmeanxjcumcor->SetBinError(i+1,hxjcumcor[i]->GetMeanError());
    hmeanxjdif->SetBinContent(i+1,hxjdif[i]->GetMean());
    hmeanxjdif->SetBinError(i+1,hxjdif[i]->GetMeanError());

    hmeanxjbkgcum->SetBinContent(i+1,hxjbkgcum[i]->GetMean());
    hmeanxjbkgcum->SetBinError(i+1,hxjbkgcum[i]->GetMeanError());
    hmeanxjbkgcumcor->SetBinContent(i+1,hxjbkgcumcor[i]->GetMean());
    hmeanxjbkgcumcor->SetBinError(i+1,hxjbkgcumcor[i]->GetMeanError());

    hmeanxjallcum->SetBinContent(i+1,hxjallcum[i]->GetMean());
    hmeanxjallcum->SetBinError(i+1,hxjallcum[i]->GetMeanError());
    hmeanxjallcumcor->SetBinContent(i+1,hxjallcumcor[i]->GetMean());
    hmeanxjallcumcor->SetBinError(i+1,hxjallcumcor[i]->GetMeanError());  

    hmeanxjeclsigallcum->SetBinContent(i+1,hxjeclsigallcum[i]->GetMean());
    hmeanxjeclsigallcum->SetBinError(i+1,hxjeclsigallcum[i]->GetMeanError());
    hmeanxjeclsigallcumcor->SetBinContent(i+1,hxjeclsigallcumcor[i]->GetMean());
    hmeanxjeclsigallcumcor->SetBinError(i+1,hxjeclsigallcumcor[i]->GetMeanError());  
    hmeanxjeclallcum->SetBinContent(i+1,hxjeclallcum[i]->GetMean());  
    hmeanxjeclallcum->SetBinError(i+1,hxjeclallcum[i]->GetMeanError());  
    hmeanxjeclallcumcor->SetBinContent(i+1,hxjeclallcumcor[i]->GetMean());  
    hmeanxjeclallcumcor->SetBinError(i+1,hxjeclallcumcor[i]->GetMeanError());    

    hmeanxjuneclsigallcum->SetBinContent(i+1,hxjuneclsigallcum[i]->GetMean());
    hmeanxjuneclsigallcum->SetBinError(i+1,hxjuneclsigallcum[i]->GetMeanError());
    hmeanxjuneclsigallcumcor->SetBinContent(i+1,hxjuneclsigallcumcor[i]->GetMean());
    hmeanxjuneclsigallcumcor->SetBinError(i+1,hxjuneclsigallcumcor[i]->GetMeanError());  
    hmeanxjuneclallcum->SetBinContent(i+1,hxjuneclallcum[i]->GetMean());  
    hmeanxjuneclallcum->SetBinError(i+1,hxjuneclallcum[i]->GetMeanError());  
    hmeanxjuneclallcumcor->SetBinContent(i+1,hxjuneclallcumcor[i]->GetMean());  
    hmeanxjuneclallcumcor->SetBinError(i+1,hxjuneclallcumcor[i]->GetMeanError());    

    hxjuneclallSUBcumcor[i]->Add(hxjuneclallcumcor[i],hxjuneclallNScumcor[i],1,-1);

    hmeanxjuneclallNScum->SetBinContent(i+1,hxjuneclallNScum[i]->GetMean());
    hmeanxjuneclallNScum->SetBinError(i+1,hxjuneclallNScum[i]->GetMeanError());
    hmeanxjuneclallNScumcor->SetBinContent(i+1,hxjuneclallNScumcor[i]->GetMean());
    hmeanxjuneclallNScumcor->SetBinError(i+1,hxjuneclallNScumcor[i]->GetMeanError());
    hmeanxjuneclallSUBcumcor->SetBinContent(i+1,hxjuneclallSUBcumcor[i]->GetMean());
    hmeanxjuneclallSUBcumcor->SetBinError(i+1,hxjuneclallSUBcumcor[i]->GetMeanError());  

  }

  vector<float > pur, combpurityAS, combpurityNS;

  for (int i=0;i<ncsv;i++) {
    pur.push_back(vhPairCode[i]->GetBinContent(1)/vhPairCode[i]->Integral());
    cout<< csvrange[i]<<" : "<<1E6*effnum[i]<<" - "<<1E6*effden<<endl;
    effnum[i]/=effden;
    combpurityAS.push_back(combpurAS[i]/combpurdenAS[i]);
    combpurityNS.push_back(combpurNS[i]/combpurdenNS[i]);
  }
  Normalize({hefficiency});

  for (auto x:pur)cout<<x<<" "; cout<<endl;

  if (MC) {
    CumDefPlot(pur,hpurity,"purity");
    CumDefPlot(effnum,hefficiency,"efficiency");
    CumDefPlot(combpurityAS,hcombpurityAS,"purity AS combinatorial");
    CumDefPlot(combpurityNS,hcombpurityNS,"purity NS combinatorial");
  }
  //purity
  // vector<TGraph *> gg;
  // auto c = getc();
  // c->SetGrid();
  // plotlegendpos = TopLeft;

  // auto l = new TLegend(0.18,0.65,0.4,0.8);
  // auto g = new TGraph(ncsv,&csvrange[0],&(pur[0]));
  // g->SetMinimum(0);
  // g->SetMaximum(1);
  // g->SetMarkerColor(kRed);//TColor::GetColorDark(i>=10 ? i+1 : i));
  // gg.push_back(g);

  // g->Draw("APL");//i==0 ? "APL" : "PL,same");
  // g->GetXaxis()->SetTitle("csv2");
  // g->GetYaxis()->SetTitle("purity");

  // hpurity->Draw("same");

  // l->SetHeader(binstr);
  // l->AddEntry(g,"cumulative","P");
  // l->AddEntry(hpurity,"differential","P");
  // l->Draw();

  // SavePlots(c,pp ? "purity_pp" : Form("purity_%d_%d",binMin,binMax));

  //efficiency
  // auto c2 = getc();
  // c2->SetGrid();
  // plotlegendpos = TopLeft;

  // auto l2 = new TLegend(0.18,0.65,0.4,0.8);
  // auto g2 = new TGraph(ncsv,&csvrange[0],&(effnum[0]));
  // g2->SetMinimum(0);
  // g2->SetMaximum(1);
  // g2->SetMarkerColor(kRed);//TColor::GetColorDark(i>=10 ? i+1 : i));

  // g2->Draw("APL");//i==0 ? "APL" : "PL,same");
  // g2->GetXaxis()->SetTitle("csv2");
  // g2->GetYaxis()->SetTitle("efficiency");

  // hefficiency->Draw("same");

  // l2->SetHeader(binstr);
  // l2->AddEntry(g2,"cumulative","P");
  // l2->AddEntry(hefficiency,"differential","P");
  // l2->Draw();

  // SavePlots(c2,pp ? "efficiency_pp" : Form("efficiency_%d_%d",binMin,binMax));


  // //combinatorial purity AS
  // auto c2 = getc();
  // c2->SetGrid();
  // plotlegendpos = TopLeft;

  // auto l2 = new TLegend(0.18,0.65,0.4,0.8);
  // auto g2 = new TGraph(ncsv,&csvrange[0],&(combpurityAS[0]));
  // g2->SetMinimum(0);
  // g2->SetMaximum(1);
  // g2->SetMarkerColor(kRed);

  // g2->Draw("APL");
  // g2->GetXaxis()->SetTitle("csv2");
  // g2->GetYaxis()->SetTitle("efficiency");

  // // hefficiency->Draw("same");

  // l2->SetHeader(binstr);
  // l2->AddEntry(g2,"cumulative","P");
  // l2->AddEntry(hefficiency,"differential","P");
  // l2->Draw();


  // //combinatorial purity NS
  // auto c2 = getc();
  // c2->SetGrid();
  // plotlegendpos = TopLeft;

  // auto l2 = new TLegend(0.18,0.65,0.4,0.8);
  // auto g2 = new TGraph(ncsv,&csvrange[0],&(combpurityNS[0]));
  // g2->SetMinimum(0);
  // g2->SetMaximum(1);
  // g2->SetMarkerColor(kRed);

  // g2->Draw("APL");
  // g2->GetXaxis()->SetTitle("csv2");
  // g2->GetYaxis()->SetTitle("efficiency");

  // // hefficiency->Draw("same");

  // l2->SetHeader(binstr);
  // l2->AddEntry(g2,"cumulative","P");
  // l2->AddEntry(hefficiency,"differential","P");
  // l2->Draw();



  plotlegendpos = None;
  aktstring = binstr;
  // Draw({hpurity});

  plotyline = hxjtrue->GetMean();
  cout<<"yline = "<<plotyline<<endl;

  plottextposx = 0.2; plottextposy = 0.81;
  plotymin = 0.64;
  plotymax = 0.71;
  plotsecondline = "cumulative";
  if (MC) Draw({hmeanxjcum});
  plotsecondline = "cumulative corrected";
  if (MC) Draw({hmeanxjcumcor});
  plotsecondline = "differential";
  if (MC) Draw({hmeanxjdif});

  plotymin = 0.5;
  plotymax = 0.8;
  plotsecondline = "background cumulative";
  if (MC) Draw({hmeanxjbkgcum});
  plotsecondline = "background cumulative corrected";
  if (MC) Draw({hmeanxjbkgcumcor});

  plotymin = 0.5;
  plotymax = 0.8;
  plotsecondline = "all cumulative";
  if (MC) Draw({hmeanxjallcum});
  plotsecondline = "all cumulative corrected";
  if (MC) Draw({hmeanxjallcumcor});


  plotyline2 = hxjeclipsedtrue->GetMean();
  cout<<"yline = "<<plotyline<<endl;


  plotsecondline = "eclipsed AS signal";
  // Draw({hmeanxjeclsigallcum});
  plotsecondline = "eclipsed AS signal corrected";
  if (MC) Draw({hmeanxjeclsigallcumcor});

  plotsecondline = "eclipsed AS all";
  Draw({hmeanxjeclallcum});
  plotsecondline = "eclipsed AS all corrected";
  Draw({hmeanxjeclallcumcor});


  plotsecondline = "UN eclipsed AS signal";
  // Draw({hmeanxjuneclsigallcum});
  plotsecondline = "UN eclipsed AS signal corrected";
  if (MC) Draw({hmeanxjuneclsigallcumcor});

  plotsecondline = "UN eclipsed AS all";
  // Draw({hmeanxjuneclallcum});
  plotsecondline = "UN eclipsed AS all corrected";
  Draw({hmeanxjuneclallcumcor});

  // plotsecondline = "UN eclipsed SUB all corrected";
  Draw({hmeanxjuneclallSUBcumcor});

  plotymax = 0.8;
  plotymin = 0;
  plotsecondline = "UN eclipsed NS all corrected";
  Draw({hmeanxjuneclallNScumcor});
  plotsecondline = "UN eclipsed NS all NOT corrected";
  Draw({hmeanxjuneclallNScum});



  //final
  if (MC && pp) plotsecondline = "Pythia 6";
  if (MC && !pp) plotsecondline = "Pythia 6 + Hydjet";
  if (!MC) plotsecondline = "Data";

  auto g = getgraph(hmeanxjuneclallSUBcumcor);

  float nomvalue=hmeanxjuneclallSUBcumcor->GetBinContent(csvrange.size()-2);
  auto c = getc();
  g->Draw("AP");//hmeanxjuneclallSUBcumcor
  g->GetXaxis()->SetRangeUser(0,1);
  g->GetXaxis()->SetTitle("Subleading jet CSV threshold");
  g->GetYaxis()->SetTitle(hmeanxjuneclallSUBcumcor->GetYaxis()->GetTitle());

  auto linenom = new TLine(0,nomvalue,1,nomvalue);
  linenom->SetLineColor(kBlack);
  //linenom->SetLineStyle(2);
  linenom->Draw();
  auto line1 = new TLine(0,nomvalue+sys,1,nomvalue+sys);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->Draw();
  auto line2 = new TLine(0,nomvalue-sys,1,nomvalue-sys);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(2);
  line2->Draw();

  cout<<"true " << hxjtrue->GetMean()<<endl;
  auto linetrue = new TLine(0,  hxjtrue->GetMean(),1,  hxjtrue->GetMean());
  linetrue->SetLineColor(kRed);
  linetrue->SetLineStyle(7);
  linetrue->Draw();



  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.2,0.8,aktstring);
  Tl->DrawLatexNDC(0.2,0.75,plotsecondline);

  SavePlot(c,TString("fullsub")+(MC ? "MC" : "Data")+(pp ? "_pp" : Form("_bin_%d_%d",binMin,binMax)));




  plotymax = 9999;
  plotymin = 9999;
  plotyline = 9999;
  plotyline2 = 9999;
}


void csvvarmc(bool MC = false)
{
  TString s = MC ? "mc" : "data";
  macro m("csvvarmc_"+s+"_0704",false);

  initcor(false);

  makeplots(-1,-1, MC,0.006);
  makeplots(60,200,MC,0.006);
  makeplots(20,60, MC,0.008);

  initcor(true);

  makeplots(0,20,  MC,0.014);
}
