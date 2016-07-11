#include "../helpers/plotting.h"
#include "../helpers/config.h"

float weightedsum(TH1F *h, int bin1, int bin2)
{
  float r = 0;
  for (int i=bin1;i<=bin2;i++)
    r+=h->GetBinCenter(i)*h->GetBinContent(i);
  return r;
}

TGraphErrors *makesysplot(TH1F *h, float dmu, int color)
{

  auto h1 = (TH1F *)h->Clone("h1");
  auto h2 = (TH1F *)h->Clone("h2");

  //if (h->GetTitle()=="Data b-jets 12")
    h->SetTitle("Data b-jets");
  h1->SetTitle("variation up");
  h2->SetTitle("variation down");

  float mu = h->GetMean();

  int n = round(mu*10);
  float i1=h->Integral(1,n);
  float w1 = weightedsum(h,1,n);

  float t = (mu-w1)/(1-i1);
  float alpha1 = (mu-dmu-t)/(w1-t*i1)-1;
  float beta1 = (1-(1+alpha1)*i1)/(1-i1)-1;
  float alpha2 = (mu+dmu-t)/(w1-t*i1)-1;
  float beta2 = (1-(1+alpha2)*i1)/(1-i1)-1;

  auto hbar = (TH1F *)h->Clone("hbar");


  for (int i=1;i<=h->GetNbinsX();i++) {
    float p = i>n ? 1+alpha1 : 1+beta1;
    h1->SetBinContent(i,h->GetBinContent(i)*p);
    p = i>n ? 1+alpha2 : 1+beta2;
    h2->SetBinContent(i,h->GetBinContent(i)*p);
  }
  Normalize({h1,h2}); //in principle it is not needed
  for (int i=1;i<=h->GetNbinsX();i++)
    hbar->SetBinError(i,fabs(h2->GetBinContent(i)-h1->GetBinContent(i))/2);



  hbar->SetFillColor(color);
  
  plotlegendpos = TopLeft;
  // plotsecondline = Form("%.3f,%.3f",h1->GetMean()-mu,mu-h2->GetMean()); //check means
  // plotthirdline = Form("%.3f",dmu);

  TGraphErrors *sys = new TGraphErrors(hbar->GetNbinsX());
  sys->SetFillColor(color);
  for (int i=1;i<=hbar->GetNbinsX();i++) {
    sys->SetPoint(i,hbar->GetBinCenter(i),hbar->GetBinContent(i));
    sys->SetPointError(i,0.05,hbar->GetBinError(i));
  }


  //draw true value, variation1 and variation2
  plotytitle = "Event fractions";
  Draw({h,h1,h2},{"E1","hist","hist","e1p"});


  // auto c2 = getc();
  // h->Draw("e1");
  // sys->Draw("e2");
  // h->Draw("e1,same");
  // SavePlot(c2,Form("%ssys2",h->GetName()));


  return sys;


  //draw cdf (for fun)
  // auto hint = h->GetIntegral();
  // auto h1int = h1->GetIntegral();
  // auto h2int = h2->GetIntegral();

  // vector<double> xaxis;
  // for (int i=1;i<=h->GetNbinsX()+1;i++) xaxis.push_back(h->GetBinLowEdge(i));

  // auto g = new TGraph(xaxis.size()+1,&xaxis[0],hint); g->SetLineColor(kRed);
  // auto g1 = new TGraph(xaxis.size()+1,&xaxis[0],h1int); g1->SetLineColor(kGreen);
  // auto g2 = new TGraph(xaxis.size()+1,&xaxis[0],h2int); g2->SetLineColor(kBlue);

  // auto c = getc();
  // g->SetMinimum(0);g->SetMaximum(1);
  // g->Draw();
  // g1->Draw("same");
  // g2->Draw("same");

  // g->GetXaxis()->SetTitle("x_{J}");
  // g->GetXaxis()->SetTitle("prob");

  // SavePlots(c,"cum");

}




void plotXJ(int inc_or_bjet=1)
{
  macro m("plotXJ_0708_default");
  
  TFile *fin = new TFile("results_0708_default/xJdphi.root");

  string species = "inc";
  if(inc_or_bjet) species = "bjt";

  TH1F *hData010 = (TH1F*) fin->Get(Form("xJ_data_%s_0_10",species.c_str()));
  TH1F *hData1030 = (TH1F*) fin->Get(Form("xJ_data_%s_10_30",species.c_str()));
  TH1F *hData30100 = (TH1F*) fin->Get(Form("xJ_data_%s_30_100",species.c_str()));
  TH1F *hDataPP = (TH1F*) fin->Get(Form("xJ_data_%s_pp",species.c_str()));
  
  TH1F *hMC010 = (TH1F*) fin->Get(Form("xJ_mc_%s_0_10",species.c_str()));
  TH1F *hMC1030 = (TH1F*) fin->Get(Form("xJ_mc_%s_10_30",species.c_str()));
  TH1F *hMC30100 = (TH1F*) fin->Get(Form("xJ_mc_%s_30_100",species.c_str()));
  TH1F *hMCPP = (TH1F*) fin->Get(Form("xJ_mc_%s_pp",species.c_str()));
  if(inc_or_bjet) hMCPP = (TH1F*) fin->Get(Form("xJ_mc_%s_pp;1",species.c_str()));

Normalize({hMC010,hMC1030,hMC30100,hMCPP}); //temporary fix

  vector<float> syserr;
  if (inc_or_bjet) syserr = {0.023,0.018,0.014,0.008}; //0-10%, 10-30%, 30-100%, pp
  else syserr = {0.023,0.016,0.010,0.007};
  int syscolor = inc_or_bjet ? kredLight : kblueLight;

  auto hData010sys =    makesysplot(hData010,syserr[0],syscolor);
  auto hData1030sys =   makesysplot(hData1030,syserr[1],syscolor);
  auto hData30100sys =  makesysplot(hData30100,syserr[2],syscolor);
  auto hDataPPsys =     makesysplot(hDataPP,syserr[3],syscolor);



  TCanvas *c1=new TCanvas("c1","c1",600,600);

  int color = kblue;
  if(inc_or_bjet) color=kred;

  
  hData010->SetMarkerColor(color);
  hData010->SetLineColor(color);
  hMC010->SetLineColor(color);

  hData010->GetXaxis()->CenterTitle(1);
  hData010->GetYaxis()->CenterTitle(1);
  hData010->SetYTitle("Event fraction");
  hData010->Draw();
  hData010sys->Draw("e2,same");
  hData010->Draw("same"); //b/c systematics should be on the back

  hMC010->SetMarkerSize(0);
  hMC010->Draw("h,same");
  
  
  CMS_lumi(c1, iPeriod, iPos ); 
    
  TLegend *l =new TLegend(0.6,0.6,0.9,0.8);
l->AddEntry(hData010,"Data","P");
l->AddEntry(hMC010,"Pythia6","l");
 if(inc_or_bjet)l->SetHeader("b-dijets");
 else l->SetHeader("Inclusive dijets");
 l->SetFillStyle(0);
 l->Draw();



  SavePlot(c1,"xJ010"+species);

  TCanvas *c2=new TCanvas("c2","c2",600,600);


  hData1030->SetMarkerColor(color);
  hData1030->SetLineColor(color);
  hMC1030->SetLineColor(color);

  hData1030->GetXaxis()->CenterTitle(1);
  hData1030->GetYaxis()->CenterTitle(1);
  hData1030->SetYTitle("Event fraction");
  hData1030->Draw();
  hData1030sys->Draw("e2,same");
  hData1030->Draw("same");
  
  hMC1030->SetMarkerSize(0);
  hMC1030->Draw("h,same");
  
  
  CMS_lumi(c2, iPeriod, iPos );

  l->Draw();



  SavePlot(c2,"xJ1030"+species);

    TCanvas *c3=new TCanvas("c3","c3",600,600);


  hData30100->SetMarkerColor(color);
  hData30100->SetLineColor(color);
  hMC30100->SetLineColor(color);

  hData30100->GetXaxis()->CenterTitle(1);
  hData30100->GetYaxis()->CenterTitle(1);
  hData30100->SetYTitle("Event fraction");
  hData30100->Draw();
  hData30100sys->Draw("e2,same");
  hData30100->Draw("same");
  
  hMC30100->SetMarkerSize(0);
  hMC30100->Draw("h,same");
  
  
  CMS_lumi(c3, iPeriod, iPos );

  l->Draw();
  
SavePlot(c3,"xJ30100"+species);



  lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp)";
  

   TCanvas *c4=new TCanvas("c4","c4",600,600);


  hDataPP->SetMarkerColor(color);
  hDataPP->SetLineColor(color);
  hMCPP->SetLineColor(color);

  hDataPP->GetXaxis()->CenterTitle(1);
  hDataPP->GetYaxis()->CenterTitle(1);
  hDataPP->SetYTitle("Event fraction");
  hDataPP->Draw();
  hDataPPsys->Draw("e2,same");  
  hDataPP->Draw("same");
  hMCPP->SetMarkerSize(0);
  hMCPP->Draw("h,same");
  
  
  CMS_lumi(c4, iPeriod, iPos );


  
  l->Draw();



SavePlot(c4,"xJpp"+species);
 
}
