#include "TF1.h"
#include "TGraph.h"

int eclipsemode = 0;
int bkgsubmode = 0;

vector<float> binbounds = {0,5,10,15,20,30,40,50,60,200};


int geteclbinindex(float bin)
{
  for(unsigned i=0;i<binbounds.size();i++)
    if (binbounds[i]>bin) return i-1;
  return binbounds.size()-2;
}


vector<float> binmeanmc, binmeandt;
vector<float> NSfracincmc, NSfracincdt;

vector<float> par0mc, par1mc;
vector<float> par0dt, par1dt;


TF1 *fecl;




TGraph *gpar0mc=0, *gpar1mc=0;
TGraph *gpar0dt=0, *gpar1dt=0;
TGraph *gNSfracincmc=0, *gNSfracincdt=0;

vector<float> constant(float v)
{
  vector<float> res;
  for (unsigned i=0;i<binbounds.size();i++) res.push_back(v);
  return res;
}

void loadeclipse()
{
  fecl = new TF1("feclipse","exp(-[0]*exp(-[1]*x))");

  if (bkgsubmode==0) {
    NSfracincmc = constant(1);
    NSfracincdt = constant(1);
  } else
  if (bkgsubmode==1) { //derived from pp. Not used!
    NSfracincmc = {0.936427,0.916281,0.899931,0.896401,0.863478,0.801669,0.684197,0.610484,0.389421};
    NSfracincdt = {0.945512,0.925795,0.910718,0.901042,0.879364,0.838917,0.771873,0.661262,0.356513};
  } else
  if (bkgsubmode==2) {
    NSfracincmc = constant(1.3);
    NSfracincdt = constant(1.3);
  } else
  if (bkgsubmode==3) { //don't subtract!
    NSfracincmc = constant(0);
    NSfracincdt = constant(0);
  }

  // //ak3
  // if (eclipsemode==0) { //correspond to bkgsubmode==0
  //   binmeanmc = {1.75025,6.93431,11.9945,17.0708,24.3864,34.4641,44.4905,54.4193,89.2678};
  //   par0mc = {19.8023,7.77764,8.33674,3.74624,1.67278,0.685874,0.903832,0.566409,0.259135};
  //   par1mc = {0.0904375,0.0808678,0.0905465,0.068902,0.0623696,0.0475709,0.0598716,0.0560104,0.0489932};

  //   binmeandt = {1.925,6.94804,12.0099,16.9759,24.3983,34.3286,44.2767,54.2765,86.1635};
  //   par0dt = {9.40264,3.88116,2.9451,1.45471,0.907856,0.445012,0.329259,0.275224,0.185311};
  //   par1dt = {0.0857002,0.0770037,0.0800895,0.0689155,0.0637999,0.0551764,0.0529313,0.0509068,0.0466953};

  // }

  // //ak4 sansjec
  // if (eclipsemode==0) { 
  //   binmeanmc = {1.84075,7.02216,11.8736,17.0287,24.0947,34.226,44.455,54.1282,92.0606};
  //   NSfracincmc = {1,1,1,1,1,1,1,1,1};
  //   par0mc = {165.371,107.106,183.464,105.364,69.4839,12.7507,4.2749,2.02672,0.549162};
  //   par1mc = {0.109613,0.106187,0.125452,0.120155,0.119234,0.0918033,0.0786657,0.0671891,0.0551698};

  //   binmeandt = {1.82588,6.8778,11.9875,16.9508,24.3868,34.3288,44.3153,54.3346,86.5521};
  //   NSfracincdt = {1,1,1,1,1,1,1,1,1};
  //   par0dt = {107.329,94.9153,79.6305,55.5985,28.8114,14.24,4.88436,1.32637,0.338164};
  //   par1dt = {0.104993,0.110982,0.112773,0.110855,0.103912,0.100325,0.0885157,0.0705909,0.0537927};
  // }
  //ak4
  if (eclipsemode==0) { //correspond to bkgsubmode==0
    binmeanmc = {1.85759,7.03006,11.9243,17.0175,24.0874,34.2351,44.4403,54.1906,92.2865};
    par0mc = {20.7652,23.7216,39.5145,24.7379,22.205,7.21922,3.56058,2.54298,0.642064};
    par1mc = {0.0799519,0.0853334,0.102825,0.0959876,0.0985177,0.0799563,0.0746687,0.0693173,0.056642};

    //with NS subid2==0
    // binmeanmc = {1.81713,7.05133,11.8931,16.9877,24.0187,34.1678,44.3589,54.1973,93.8959};
    // par0mc = {23.3885,28.7707,55.9523,52.7231,44.0016,10.3153,10.9794,6.14483,1.46164};
    // par1mc = {0.0832596,0.0914355,0.113106,0.115699,0.117081,0.0938252,0.111851,0.104616,0.105466};

    
    binmeandt = {1.84443,6.8932,11.9979,16.9695,24.4043,34.3534,44.3177,54.3508,86.6278};
    par0dt = {19.5221,18.4325,21.1487,17.3539,12.0889,8.20492,3.8046,1.342,0.371563};
    par1dt = {0.08256,0.0877552,0.0934504,0.0928299,0.0893618,0.0895741,0.0824222,0.069007,0.0544498};

  }


  if (eclipsemode==1) { //correspond to bkgsubmode==1

    binmeanmc = {1.85759,7.03006,11.9243,17.0175,24.0874,34.2351,44.4403,54.1906,92.2865};
    NSfracincmc = {0.936427,0.916281,0.899931,0.896401,0.863478,0.801669,0.684197,0.610484,0.389421};
    par0mc = {18.1711,19.9312,32.1561,20.4474,17.6702,5.41424,2.3534,1.46491,0.244543};
    par1mc = {0.0789364,0.0839883,0.10115,0.0946377,0.0971702,0.0789454,0.0742867,0.068507,0.0566169};


    binmeandt = {1.84443,6.8932,11.9979,16.9695,24.4043,34.3534,44.3177,54.3508,86.6278};
    NSfracincdt = {0.945512,0.925795,0.910718,0.901042,0.879364,0.838917,0.771873,0.661262,0.356513};
    par0dt = {17.5464,15.8957,17.9866,14.7832,10.095,6.53038,2.84406,0.874138,0.14021};
    par1dt = {0.0817867,0.0864993,0.0923211,0.091923,0.0885183,0.0887129,0.0820349,0.0689021,0.0556841};

  }
  // if (eclipsemode==2) {
  //   binmeanmc = {1.83014,7.02054,11.8852,17.0213,24.1649,34.2033,44.4239,54.1938,92.1712};
  //   NSfracincmc = {1.04845,1.07879,1.09039,1.097,1.10314,1.18029,1.34007,1.42532,1.51716};
  //   par0mc = {45.4145,38.8973,49.1055,26.7811,45.9855,11.6104,4.98896,2.9713,1.12204};
  //   par1mc = {0.0898231,0.0910681,0.105007,0.0983384,0.108087,0.083725,0.0754124,0.0666132,0.0584176};

  //   binmeandt = {1.82676,6.8791,11.9878,16.9495,24.465,34.3395,44.3056,54.3455,86.5526};
  //   NSfracincdt = {1.04607,1.06768,1.09133,1.11049,1.10314,1.14909,1.23484,1.37751,1.61095};
  //   par0dt = {36.9319,28.4173,26.2258,17.1042,20.3553,11.2449,3.97173,1.35135,0.647808};
  //   par1dt = {0.0895161,0.0925859,0.0959309,0.0921781,0.0947111,0.09177,0.079296,0.0641041,0.055368};


  // }
  // if (eclipsemode==3) {
  //   binmeanmc = {1.83014,7.02054,11.8852,17.0213,24.1649,34.2033,44.4239,54.1938,92.1712};
  //   NSfracincmc = {0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};
  //   par0mc = {18.0721,15.9602,21.6177,13.351,20.9457,5.62291,2.31539,1.3144,0.479119};
  //   par1mc = {0.0818099,0.0838429,0.098696,0.0943111,0.102402,0.0805474,0.0736999,0.0651855,0.0572982};


  //   binmeandt = {1.82676,6.8791,11.9878,16.9495,24.465,34.3395,44.3056,54.3455,86.5526};
  //   NSfracincdt = {0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};
  //   par0dt = {16.0205,13.612,12.8519,8.78273,10.4179,5.80439,2.07874,0.645129,0.284836};
  //   par1dt = {0.0827475,0.0875804,0.0915125,0.0888582,0.0912506,0.0889744,0.0780994,0.0632023,0.0558225};


  // }
  // if (eclipsemode==4) {
  //   binmeanmc = {1.83014,7.02054,11.8852,17.0213,24.1649,34.2033,44.4239,54.1938,92.1712};
  //   NSfracincmc = {1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3};
  //   par0mc = {90.9637,65.7384,74.0837,36.6212,64.5962,13.4575,4.81417,2.63853,0.956353};
  //   par1mc = {0.0973245,0.0964959,0.108977,0.100716,0.111066,0.0845123,0.0753898,0.0661916,0.058392};

  //   binmeandt = {1.82676,6.8791,11.9878,16.9495,24.465,34.3395,44.3056,54.3455,86.5526};
  //   NSfracincdt = {1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3};
  //   par0dt = {69.0783,43.7725,36.6295,22.0801,27.5788,13.3352,4.21767,1.27131,0.515417};
  //   par1dt = {0.0960085,0.0963982,0.0985108,0.0937456,0.0970974,0.0925128,0.0794284,0.064097,0.0551593};

  // }



  gpar0mc = new TGraph(binmeanmc.size(),&binmeanmc[0],&par0mc[0]);
  gpar1mc = new TGraph(binmeanmc.size(),&binmeanmc[0],&par1mc[0]);

  gpar0dt = new TGraph(binmeandt.size(),&binmeandt[0],&par0dt[0]);
  gpar1dt = new TGraph(binmeandt.size(),&binmeandt[0],&par1dt[0]);

  gNSfracincmc = new TGraph(binmeanmc.size(),&binmeanmc[0],&NSfracincmc[0]);
  gNSfracincdt = new TGraph(binmeandt.size(),&binmeandt[0],&NSfracincdt[0]);

}

float eclipseWeightmc(float jtpt2, float bin)
{
  if (gpar0mc==0)  loadeclipse();

  if (jtpt2<40) return 0;

  float par0 = gpar0mc->Eval(bin);
  float par1 = gpar1mc->Eval(bin);


  fecl->SetParameters(par0,par1);
  float p = fecl->Eval(jtpt2);

  return 1/p;
}

float eclipseWeightdt(float jtpt2, float bin)
{
  if (gpar0dt==0)  loadeclipse();
  int b = geteclbinindex(bin);

  float par0 = gpar0dt->Eval(bin);
  float par1 = gpar1dt->Eval(bin);

  fecl->SetParameters(par0,par1);
  float p = fecl->Eval(jtpt2);

  return 1/p;
}

float NSfracmc(float bin)
{
  if (gNSfracincmc==0)  loadeclipse();
  return gNSfracincmc->Eval(bin);
}

float NSfracdt(float bin)
{
  if (gNSfracincdt==0)  loadeclipse();
  return gNSfracincdt->Eval(bin);
}

void drawNSfractions()
{
  eclipsemode = 1;
  loadeclipse();

  gNSfracincmc->SetMarkerColor(kBlue);
  gNSfracincdt->SetMarkerColor(kRed);
  gNSfracincmc->SetMinimum(0);
  gNSfracincmc->SetMaximum(1);
  auto c = new TCanvas("c","c",600,600);
  gNSfracincmc->Draw("AP");
  gNSfracincdt->Draw("P,same");
  auto l = new TLegend(0.2,0.2,0.4,0.4);
  l->AddEntry(gNSfracincmc,"Pythia 6 + Hydjet","P");
  l->AddEntry(gNSfracincdt,"Data","P");
  l->Draw();
  gNSfracincmc->GetXaxis()->SetTitle("bin");
  gNSfracincmc->GetYaxis()->SetTitle("combinatorial background fraction");

  c->SaveAs("NSfractions.pdf");
}