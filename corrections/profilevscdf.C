#include "../helpers/plotting.h"
#include "../helpers/config.h"
#include "TProfile.h"
// #include "eclipsederive.C"

vector<int> binbounds = {0,5,10,15,20,30,40,50,60,200};

TGraph *getGraph(TH1 *p)
{
  vector<double>vx,vy;
  for (int i=1;i<=p->GetNbinsX();i++) {
    vx.push_back(p->GetBinCenter(i));
    vy.push_back(p->GetBinContent(i));
  }
  return new TGraph(p->GetNbinsX(), &vx[0],&vy[0]);
}

float findx(TGraph *g, float s)
{
  float xmin = 0, xmax = 300;
  float r =(xmax+xmin)/2; 
  float y=g->Eval(r);
  float eps = 0.01;
  while (abs(y-s)>eps && abs(xmax-xmin)>eps) {
    float d = xmax - xmin;
    if (g->Eval(xmax-d/2)>s) xmax-=d/2;
    if (g->Eval(xmin+d/2)<s) xmin+=d/2;
    r=(xmax+xmin)/2; 
    y=g->Eval(r);
  }
  return r;
}

double getMedian(TH1 *h)
{
  double m[1];
  double p[1] = {0.5};
  h->GetQuantiles(1,m,p);
  return m[0];
}

double histdist(TH1F *h1, TH1F *h2)
{
  double s = 0;
  int n=h1->GetNbinsX();
  for (int i=1;i<=n;i++) {
    double t = h1->GetBinContent(i)-h2->GetBinContent(i);
    s+=t*t;
  }

  return s;
}

TH1F *shifthist(TH1F *h, int nb) 
{
  auto ht=(TH1F *)h->Clone(Form("ht%d",nb));
  for (int i=1;i<=h->GetNbinsX();i++)
    ht->SetBinContent(i,h->GetBinContent(i+nb));
  return ht;
}

float finddelta(TH1F *h1, TH1F *h2)
{
  Normalize({h1,h2});
  int d = -10;
  int dmin = -10;
  double distmin = 99999;
  while (d<10) {
    d++;
    auto shiftedh = shifthist(h2,d);
    double t = histdist(h1,shiftedh);
    if (t<distmin) {
      distmin = t;
      dmin = d;
    }
    Draw({h1,shiftedh});
    cout<<d<<" "<<t<<endl;
  }
  return dmin*h1->GetBinWidth(1);
}

void ScaleVisibleBins(TH1F *h, float SF)
{
  float underv = h->GetBinContent(0);
  float undere = h->GetBinError(0);

  float i = h->Integral();
  h->Scale(SF);

  h->SetBinContent(0,underv+i*(1-SF));
  h->SetBinError(0,undere+i*(1-SF));
}

vector<TH1F *> eclipsecurves;
vector<TString> eclipsecurvenames;

float runbins(int b1, int b2) 
{

  auto fdt = config.getfile_djt("dtPbjcl");
  auto ntdt = (TTree *)fdt->Get("nt");

  seth(71,38,180);
  auto hdt = geth(Form("hdt%d%d",b1,b2),"Data; Max jet p_{T} [GeV]; Event fractions");
 
  ntdt->Project(hdt->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05)",b1,b2));
  hdt->SetBinContent(1,hdt->GetBinContent(0)+hdt->GetBinContent(1));
  auto gdt = getCDF(hdt);
  SetData({hdt});
  // Draw({gdt});
  // for (int i=0;i<gdt->GetNbinsX();i++)
    // gdt->SetBinError(i,0);


  auto file = config.getfile_djt("mcPbqcd");
  auto nt = (TTree *)file->Get("nt");

  int pthatcut = 80;


  auto p = new TProfile(Form("p%d%d",b1,b2),"prof;p_{T,2} threshold [GeV];subleading jet findiding efficiency",35,40,180);

  
  nt->Project(p->GetName(),"(subid2 == 0 && refpt2 > 20):jtptSignal2",Form("weight*(jtpt1>100 &&bin>=%d && bin<%d && pthat>%d)",b1,b2,pthatcut));


  auto h = geth(Form("h%d%d",b1,b2),"h all");
  auto h2 = geth(Form("h2%d%d",b1,b2),"h not signal");
  nt->Project(h->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05 && pthat>%d)",b1,b2,pthatcut));
  //to check the signal effect on the result
  nt->Project(h2->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d &&dphi21<1.05 && pthat>%d && !(subid2 == 0 && refpt2 > 20))",b1,b2,pthatcut));

  float SF = h2->Integral()/h->Integral();
  cout<<"h2/h integral : "<<SF<<endl;

  ScaleVisibleBins(h,SF);
  // Normalize({h,h2});
  Draw({h2,h});


  h->SetBinContent(1,(h->GetBinContent(0)+h->GetBinContent(1)));
  h2->SetBinContent(1,h2->GetBinContent(0)+h2->GetBinContent(1));

  cout<<"h2/h underbin : "<<h2->GetBinContent(1)/h->GetBinContent(1)<<endl;

  // cout<<" bins to shift "<<finddelta(h,h2)<<endl;

  auto g = getCDF(h);
  g->SetLineColor(darkgreen);
  g->SetLineWidth(3);

  eclipsecurves.push_back(g);
  eclipsecurvenames.push_back(FloatToStr(b1/2.)+"-"+FloatToStr(b2/2.)+"%");

  auto g2 = getCDF(h2);
  g2->SetLineColor(darkblue);
  g2->SetLineWidth(3);

  gdt->SetLineColor(darkred);
  gdt->SetLineWidth(3); 

  float profilemedian = findx(getGraph(p),0.5);
  float cdfmedian = findx(getCDFgraph(g),0.5);
  float cdfbkgmedian = findx(getCDFgraph(g2),0.5);//getMedian(h);

  cout<<"medians"<<endl;
  cout<<"profile   : "<<profilemedian<<endl;
  cout<<"cdf NS    : "<<getMedian(h)<<endl;
  cout<<"cdf NS bkg: "<<getMedian(h2)<<endl;

  p->SetMinimum(0);
  p->SetMaximum(1);

  auto c = getc();
  p->Draw();
  g->Draw("hist,same");
  // g2->Draw("hist,same");
  gdt->Draw("hist,same");

  p->GetXaxis()->CenterTitle();
  p->GetYaxis()->CenterTitle();


  TString centstr = FloatToStr(b1/2.)+"-"+FloatToStr(b2/2.)+"%";

  plotlegendpos = BottomRight;
  // plotlegenddx = -0.1;
  auto l=getLegend();
  l->SetY1(0.3); l->SetY2(0.5);
  l->SetHeader("Pythia 6 + Hydjet ");
  l->AddEntry(p,"true","P");//subleading jet efficiency
  // l->AddEntry(g2,"c.d.f. of bkg only jets","L");
  l->AddEntry(g,"estimated","L");//c.d.f. of near-side jets

  auto l2 = getLegend();
  l2->SetHeader("PbPb data ");
  l2->SetY1(0.2); l2->SetY2(0.3);
  l2->AddEntry(gdt,"estimated","L");//c.d.f. of near-side jets

  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.49, 0.52 , centstr);

  l->Draw();
  l2->Draw();


  CMS_lumi(c, iPeriod, 33 ); 

  SavePlot(c,Form("profilevscdf%d%d",b1,b2));

  return cdfbkgmedian-cdfmedian;//profilemedian-cdfmedian;
}

//DO NOT USE THIS ONE, but from eclipseclosure.C
void draweclipsecurves()
{
  auto c = getc();
  plotlegendpos = BottomRight;
  auto l = new TLegend(0.5,0.3,0.8,0.7);
  for (int i=0;i<eclipsecurves.size();i++) {
    eclipsecurves[i]->SetLineColorAlpha(kgreen,1-0.9*((float)i)/eclipsecurves.size());
    l->AddEntry(eclipsecurves[i],eclipsecurvenames[i],"L");
    eclipsecurves[i]->SetMaximum(1);
    eclipsecurves[i]->Draw(i==0 ? "" : "same");

    eclipsecurves[i]->GetXaxis()->SetTitle("p_{T,2} threshold [GeV]");
    eclipsecurves[i]->GetYaxis()->SetTitle("found fraction");
  }
  l->Draw();

  SavePlot(c,"eclipsecurves");
}

void profilevscdf()
{  
  macro m("profilevscdf_0704");

  // binbounds = {0,20,60,200};//
  binbounds = {0,5,10,15,20,30,40,50,60,200};

  vector<float> diff;
  for(unsigned i=0;i<binbounds.size()-1;i++) 
    diff.push_back(runbins(binbounds[i],binbounds[i+1]));

//output
  for(unsigned i=0;i<binbounds.size()-1;i++) 
    cout<<binbounds[i]/2<<" - "<<binbounds[i+1]/2<<"% : "<<diff[i]<<endl;


  draweclipsecurves();
}



