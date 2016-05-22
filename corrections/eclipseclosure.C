#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "TProfile.h"


//from profile
// vector<float> binbounds = {0,5,10,15,20,30,40,50,60,200};
// vector<float> par0mc = {352.239,560.346,392.329,343.754,256.353,243.9,146.921,133.67,40.5061};
// vector<float> par1mc = {0.0974531,0.116765,0.121361,0.127222,0.132043,0.148747,0.149254,0.157618,0.153144};
// vector<float> par0mc = {541.981,692.366,471.802,313.334,252.521,302.789,160.788,146.241,34.5569};
// vector<float> par1mc = {0.102214,0.11914,0.12276,0.123191,0.13009,0.153598,0.150147,0.162346,0.1488};
//from NS
// vector<float> binbounds = {0,5,10,15,20,30,40,50,60,200};
// vector<float> par0 = {242.885,143.471,259.942,166.572,130.455,73.0904,64.5492,48.8908,23.4887};
// vector<float> par1 = {0.0891711,0.0888366,0.108133,0.109136,0.114936,0.11348,0.121985,0.127472,0.121776};


vector<float> binbounds = {0,5,10,15,20,30,40,50,60,200};
//NS MC pthat>65
// vector<float> par0mc = {242.708,143.103,259.991,153.735,114.992,48.0445,20.631,4.39078,1.78779};
// vector<float> par1mc = {0.0891591,0.0888055,0.108135,0.107633,0.112421,0.104633,0.097206,0.0735431,0.0563029};

//NS MC pthat>50
vector<float> binmean = {1.76756, 7.02132,  11.8198, 17.0293, 23.8897, 34.1438, 44.1978, 54.1609,  91.5871};

vector<float> par0mc = {310.123,  257.007,  338.152, 278.773, 156.73,  39.3968, 30.3697, 7.18124,  2.24178};
vector<float> par1mc = {0.0931636,0.0994282,0.114211,0.120011,0.118508,0.101331,0.107097,0.0830526,0.0612657};
// NS data  
vector<float> par0dt = {271.692,254.36,196.266,187.834,107.866,67.3588,42.9682,16.8285,5.79053};
vector<float> par1dt = {0.0949912,0.104842,0.107888,0.114252,0.111524,0.109608,0.107457,0.0888285,0.0624925};

auto fecl = new TF1("feclipse","exp(-[0]*exp(-[1]*x))");

int geteclbinindex(float bin)
{
  for(unsigned i=0;i<binbounds.size();i++)
    if (binbounds[i]>bin) return i-1;
  return binbounds.size()-2;
}


float eclipseWeightdt(float jtpt2, float bin)
{
  // return 1;
  int b = geteclbinindex(bin);
  fecl->SetParameters(par0dt[b],par1dt[b]);
  float p = fecl->Eval(jtpt2);
  if (p<0.05) p=0.05;
  return 1/p;
}


TGraph *gpar0mc=0, *gpar1mc=0;
void loadcorr()
{

  //for the moment, binbounds work better than binmeans
  gpar0mc = new TGraph(binbounds.size(),&binbounds[0],&par0mc[0]);
  gpar1mc = new TGraph(binbounds.size(),&binbounds[0],&par1mc[0]);
}

float eclipseWeightmc(float jtpt2, float bin)
{
  if (gpar0mc==0)  loadcorr();

  // return 1;
  int b = geteclbinindex(bin);

// float par0 = par0mc[b];
// float par1 = par1mc[b];

  float par0 = gpar0mc->Eval(bin);
  float par1 = gpar1mc->Eval(bin);


  fecl->SetParameters(par0,par1);
  float p = fecl->Eval(jtpt2);
  if (p<0.05) p=0.05;
  return 1/p;
}

float eclipseWeightmcNoCutoff(float jtpt2, float bin)
{
  // return 1;
  int b = geteclbinindex(bin);
  // float par0 = par0mc[b];
  // float par1 = par1mc[b];

  float par0 = gpar0mc->Eval(bin);
  float par1 = gpar1mc->Eval(bin);
  fecl->SetParameters(par0,par1);
  float p = fecl->Eval(jtpt2);
  return 1/p;
}

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
  SavePlots(c,filename);

}


void checkclosure(int binMin, int binMax)
{

  auto f = config.getfile_djt("mcPbqcd");

  seth(10,0,1);
  auto hxjASsiguncorr = geth("hxjASsiguncorr","Signal, eclipsed;x_{J};Event fractions");
  auto hxjASsigcorr = geth("hxjASsigcorr","Signal, corrected;x_{J};Event fractions");

  auto hxjAStrue = geth("hxjAStrue","Signal, not eclipsed;x_{J};Event fractions");
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

  Fill(f,{"pthat","weight","jtpt2","bin","jtpt1","jtpt2","subid2","refpt2","dphi21","jtptSignal2","dphi21","dphiSignal21"}, [&](dict &d) {
    if (d["pthat"]<pthatcut) return;
    float w = d["weight"];
    float corr = eclipseWeightmc(d["jtpt2"],d["bin"]);
    hcorr->Fill(eclipseWeightmcNoCutoff(d["jtpt2"],d["bin"]),w);
    hcorr2d->Fill(eclipseWeightmcNoCutoff(d["jtpt2"],d["bin"]),d["bin"],w);
    if (d["bin"]<binMin || d["bin"]>=binMax) return;

    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut) {
      if (d["subid2"]==0 && d["refpt2"]>20 && d["dphi21"]>PI23) 
      {
        hxjASsiguncorr->Fill(d["jtpt2"]/d["jtpt1"],w);
        hxjASsigcorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr);
      }          

      if (d["dphi21"]>PI23) {
        hxjASuncorr->Fill(d["jtpt2"]/d["jtpt1"],w);
        hxjAScorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr);
      }

      if (d["dphi21"]<PI13){
        hxjNSuncorr->Fill(d["jtpt2"]/d["jtpt1"],w);
        hxjNScorr->Fill(d["jtpt2"]/d["jtpt1"],w*corr);
      }
    }

    if (d["jtpt1"]>pt1cut && d["jtptSignal2"]>pt2cut &&  d["dphiSignal21"]>PI23)
      hxjAStrue->Fill(d["jtptSignal2"]/d["jtpt1"],w);
  });
 
hxjASclos->Add(hxjNScorr,hxjAScorr,1,-1);
hxjbkgsubuncorr->Add(hxjNSuncorr,hxjASuncorr,1,-1);

NormalizeAllHists();

plotylog = true;
Draw({hcorr});
plotylog = false;

auto c = getc();
hcorr2d->Draw("colz");
c->SetLogz();
c->SetLogx();
SavePlots(c,"hcorr2d");


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


}

void eclipseclosure()
{
  macro m("eclipseclosureNS_0519");


  checkclosure(0,20);
  checkclosure(20,60);
  checkclosure(60,200);

}