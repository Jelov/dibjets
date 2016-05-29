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

// wrong - without underflow
vector<float> par0mc = {310.123,  257.007,  338.152, 278.773, 156.73,  39.3968, 30.3697, 7.18124,  2.24178};
vector<float> par1mc = {0.0931636,0.0994282,0.114211,0.120011,0.118508,0.101331,0.107097,0.0830526,0.0612657};

//???correct - with underflow bin
// //vector<float> par0mc = {340.855,285.031,379.29,312.579,176.524,43.3339,33.695,7.64328,1.10267};
// //vector<float> par1mc = {0.0931846,0.0994929,0.114226,0.1199,0.118522,0.101184,0.107032,0.0827885,0.059788};


//UPDATED WITH UNDERFLOW BIN FIX!
// vector<float> par0mc = {282.615,233.15,303.303,248.589,135.641,26.4813,19.1192,5.8134,0.89925};
// vector<float> par1mc = {0.0931685,0.0994672,0.114318,0.120107,0.117998,0.0954307,0.0998593,0.0806762,0.0583244};

//plus 1 GeV
// vector<float> par0mc = {307.678,260.455,345.329,293.864,153.239,35.5853,24.7877,6.85737,0.984627};
// vector<float> par1mc = {0.0930313,0.0997206,0.114631,0.121013,0.118095,0.0992568,0.10326,0.0818517,0.0587724};

//minus 1 GeV
// vector<float> par0mc = {255.119,213.95,277.18,225.139,122.169,19.9603,12.4845,5.3611,0.791768};
// vector<float> par1mc = {0.0930073,0.099771,0.114794,0.120535,0.118287,0.0917484,0.0935159,0.0802547,0.0569541};



// NS data -- old, without underflow
vector<float> par0dt = {271.692,254.36,196.266,187.834,107.866,67.3588,42.9682,16.8285,5.79053};
vector<float> par1dt = {0.0949912,0.104842,0.107888,0.114252,0.111524,0.109608,0.107457,0.0888285,0.0624925};

//Data - UPDATED!!!
// vector<float> par0dt = {247.354,229.656,175.195,160.035,70.5693,29.6137,10.4073,3.33657,0.417425};
// vector<float> par1dt = {0.0950106,0.104889,0.107828,0.113688,0.107395,0.101935,0.0942873,0.08042,0.054413};

//min 1 GeV
// vector<float> par0dt = {224.496,205.743,157.371,147.157,60.5557,26.7927,7.98273,2.9463,0.368765};
// vector<float> par1dt = {0.0950121,0.104765,0.107892,0.114274,0.106524,0.102094,0.0908525,0.0795072,0.0532526};

//plus 1 GeV
// vector<float> par0dt = {271.038,251.323,195.291,180.785,79.2452,37.3225,14.0793,3.87185,0.479263};
// vector<float> par1dt = {0.0949957,0.104601,0.107896,0.113841,0.107592,0.104523,0.0983239,0.0816622,0.0559585};


auto fecl = new TF1("feclipse","exp(-[0]*exp(-[1]*x))");

int geteclbinindex(float bin)
{
  for(unsigned i=0;i<binbounds.size();i++)
    if (binbounds[i]>bin) return i-1;
  return binbounds.size()-2;
}


/////OLD CODE!!!!!
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
/////OLD CODE!!!!!

//NEW CODE!!!!!
// TGraph *gpar0mc=0, *gpar1mc=0;
// TGraph *gpar0dt=0, *gpar1dt=0;
// void loadcorr()
// {

//   //for the moment, binbounds work better than binmeans
//   gpar0mc = new TGraph(binmean.size(),&binmean[0],&par0mc[0]);//binbounds
//   gpar1mc = new TGraph(binmean.size(),&binmean[0],&par1mc[0]);//binmean

//   gpar0dt = new TGraph(binmean.size(),&binmean[0],&par0dt[0]);
//   gpar1dt = new TGraph(binmean.size(),&binmean[0],&par1dt[0]);

// }

// float eclipseWeightmc(float jtpt2, float bin)
// {
//   if (gpar0mc==0)  loadcorr();

//   // return 1;
//   int b = geteclbinindex(bin);

// // float par0 = par0mc[b];
// // float par1 = par1mc[b];

//   float par0 = gpar0mc->Eval(bin);
//   float par1 = gpar1mc->Eval(bin);


//   fecl->SetParameters(par0,par1);
//   float p = fecl->Eval(jtpt2);
//   if (p<0.05) p=0.05;
//   return 1/p;
// }


// float eclipseWeightdt(float jtpt2, float bin)
// {
//   if (gpar0dt==0)  loadcorr();
//   int b = geteclbinindex(bin);

//   float par0 = gpar0dt->Eval(bin);
//   float par1 = gpar1dt->Eval(bin);

//   fecl->SetParameters(par0,par1);
//   float p = fecl->Eval(jtpt2);
//   if (p<0.05) p=0.05;
//   return 1/p;
// }


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

vector<float> xjtrue, xjsubuncorr,xjunsub,xjsubcorr;
vector<float> exjtrue, exjsubuncorr,exjunsub,exjsubcorr;

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


xjtrue.push_back(hxjAStrue->GetMean());
exjtrue.push_back(hxjAStrue->GetMeanError());

xjsubuncorr.push_back(hxjbkgsubuncorr->GetMean());
exjsubuncorr.push_back(hxjbkgsubuncorr->GetMeanError());

xjunsub.push_back(hxjASuncorr->GetMean());
exjunsub.push_back(hxjASuncorr->GetMeanError());

xjsubcorr.push_back(hxjASclos->GetMean());
exjsubcorr.push_back(hxjASclos->GetMeanError());


}

void fillhistbins(TH1F *h, vector<float> x, vector<float> e)
{
  for (unsigned i=0;i<x.size();i++) {
    h->SetBinContent(i+1,x[i]);
    h->SetBinError(i+1,e[i]);
  }
}

void eclipseclosure()
{
  macro m("eclipseclosure");


  checkclosure(0,20);
  checkclosure(20,60);
  checkclosure(60,200);

  seth(3,0,3);
  auto hxjtrue = geth("hxjtrue","Truth");
  auto hxjunsub = geth("hxjunsub","Raw");
  auto hxjsubuncorr = geth("hxjsubuncorr","Sideband subtracted");
  auto hxjsubcorr = geth("hxjsubcorr","+ Eclipse correction");

  fillhistbins(hxjtrue, xjtrue, exjtrue);
  fillhistbins(hxjsubuncorr, xjsubuncorr, exjsubuncorr);
  fillhistbins(hxjunsub, xjunsub, exjunsub);
  fillhistbins(hxjsubcorr, xjsubcorr, exjsubcorr);


  vector<TString> axnames;
 for (auto s:binnames) axnames.push_back(s);
  RenameBinLabelsX(hxjtrue,axnames);
  RenameBinLabelsX(hxjsubuncorr,axnames);
  RenameBinLabelsX(hxjunsub,axnames);
  RenameBinLabelsX(hxjsubcorr,axnames);  


  plotlegendpos = TopLeft;
  // Draw(hmcPbdphi);

  aktstring = "";

  ShuffleBins(hxjtrue,{3,2,1});
  ShuffleBins(hxjsubuncorr,{3,2,1});
  ShuffleBins(hxjunsub,{3,2,1});
  ShuffleBins(hxjsubcorr,{3,2,1});

  plotymin = 0.5;
  plotymax = 0.8;

  plotputmean = false;

SetInc({hxjtrue,hxjsubuncorr,hxjunsub,hxjsubcorr});
hxjtrue->SetMarkerColor(kBlack);
hxjsubuncorr->SetMarkerColor(kBlack);

hxjtrue->SetMarkerSize(2.);
hxjsubuncorr->SetMarkerSize(1.4);
hxjunsub->SetMarkerSize(1.4);
hxjsubcorr->SetMarkerSize(1.4);


hxjtrue->SetMarkerStyle(kOpenDiamond);
hxjsubuncorr->SetMarkerStyle(kOpenSquare);
hxjunsub->SetMarkerStyle(kOpenCircle);
hxjsubcorr->SetMarkerStyle(kFullCircle);

  hxjsubcorr->SetMinimum(0.55);
  hxjsubcorr->SetMaximum(0.8);

  auto c=getc();
  auto l = getLegend();
  l->AddEntry(hxjtrue,hxjtrue->GetTitle(),"P");
  l->AddEntry(hxjunsub,hxjunsub->GetTitle(),"P");
  l->AddEntry(hxjsubuncorr,hxjsubuncorr->GetTitle(),"P");
  l->AddEntry(hxjsubcorr,hxjsubcorr->GetTitle(),"P");


  hxjsubcorr->Draw();
  hxjtrue->Draw("same");
  hxjsubuncorr->Draw("same");

  hxjunsub->Draw("same");
  l->Draw();


  SavePlots(c,"eclipsedemonstration");
  // Draw({hxjtrue,hxjsubuncorr,hxjunsub,hxjsubcorr});

}