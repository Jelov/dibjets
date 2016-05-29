#include "../helpers/plotting.h"
#include "../helpers/config.h"
#include "TProfile.h"
#include "eclipsederive.C"


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

float runbins(int b1, int b2) 
{

  auto fdt = config.getfile_djt("dtPbj60");
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

  int pthatcut = 50;


  auto p = new TProfile(Form("p%d%d",b1,b2),"prof;p_{T,2} threshold [GeV];found fraction",70,40,180);

  
  nt->Project(p->GetName(),"(subid2 == 0 && refpt2 > 20):jtptSignal2",Form("weight*(jtpt1>100 &&bin>=%d && bin<%d && pthat>%d)",b1,b2,pthatcut));


  auto h = geth(Form("h%d%d",b1,b2),"h");
  auto h2 = geth(Form("h2%d%d",b1,b2),"h");
  nt->Project(h->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05 && pthat>%d)",b1,b2,pthatcut));
  //to check that signal has no effect on the result
  nt->Project(h2->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d &&dphi21<1.05 && pthat>%d && !(subid2 == 0 && refpt2 > 20))",b1,b2,pthatcut));

  h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
  h2->SetBinContent(1,h2->GetBinContent(0)+h2->GetBinContent(1));


  auto g = getCDF(h);
  g->SetLineColor(darkgreen);
  g->SetLineWidth(3);

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
  cout<<"cdf NS    : "<<cdfmedian<<endl;
  cout<<"cdf NS bkg: "<<getMedian(h2)<<endl;

  p->SetMinimum(0);
  p->SetMaximum(1);

  auto c = getc();
  p->Draw();
  g->Draw("hist,same");
  // g2->Draw("hist,same");
  gdt->Draw("hist,same");
  plotlegendpos = BottomRight;
  plotlegenddx = -0.1;
  auto l=getLegend();
  l->SetY1(0.3); l->SetY2(0.5);
  l->SetHeader("MC");
  l->AddEntry(p,"subleading jet efficiency","P");
  // l->AddEntry(g2,"c.d.f. of bkg only jets","L");
  l->AddEntry(g,"c.d.f. of near-side jets","L");

  auto l2 = getLegend();
  l2->SetHeader("Data");
  l2->SetY1(0.2); l2->SetY2(0.3);
  l2->AddEntry(gdt,"c.d.f. of near-side jets","L");

  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.55, 0.8, "PbPb "+FloatToStr(b1/2.)+"-"+FloatToStr(b2/2.)+"%");

  l->Draw();
  l2->Draw();



  SavePlots(c,Form("profilevscdf%d%d",b1,b2));

  return cdfbkgmedian-cdfmedian;//profilemedian-cdfmedian;
}

void profilevscdf()
{  
  macro m("profilevscdf");

  vector<float> diff;
  for(int i=0;i<binbounds.size()-1;i++) 
    diff.push_back(runbins(binbounds[i],binbounds[i+1]));

//output
  for(int i=0;i<binbounds.size()-1;i++) 
    cout<<binbounds[i]/2<<" - "<<binbounds[i+1]/2<<"% : "<<diff[i]<<endl;

}



