#ifndef PLOTTINH_H
#define PLOTTINH_H

#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLine.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TSystem.h"

using namespace std;

enum LegendPosition {TopRight, TopLeft, BottomRight, BottomLeft, None};

int darkred = TColor::GetColorDark(2);
int darkgreen = TColor::GetColorDark(3);
int darkblue = TColor::GetColorDark(4);
int lightblue = TColor::GetColorBright(4);

bool buildFromVector = false;
vector<float> buildvector;
int buildnbins = 20;
float buildxmin = 0;
float buildxmax = 1;
int buildnbinsy = 20;
float buildymin = 0;
float buildymax = 1;
vector<TH1 *> allhists;


int ccounter = 0;
TCanvas *getc()
{
  ccounter++;
  TString name = TString("c")+TString::Itoa(ccounter,10);

  return new TCanvas(name,name,600,600);

}

void buildh(int nbins=20, float xmin=0, float xmax=1, int nbinsy=20, float ymin=0, float ymax=1)
{
  buildFromVector = false;
  buildnbins = nbins;
  buildxmin = xmin;
  buildxmax = xmax;
  buildnbinsy = nbinsy;
  buildymin = ymin;
  buildymax = ymax;


}

bool firstRunMacro = true;
TString histoutputfilename;

void buildh(vector<float> &binsvector)
{
  buildFromVector = true;
  buildvector = binsvector;
}

TString buildNamesuffix = "";
TString buildTitlesuffix = "";


TH1F *geth(TString hname, TString htitle)
{
  if (firstRunMacro) {
    TH1F *h;
    if (buildFromVector)
      h = new TH1F(hname+buildNamesuffix, buildTitlesuffix+htitle,buildvector.size()-1,&buildvector[0]);
    else
      h = new TH1F(hname+buildNamesuffix, buildTitlesuffix+htitle, buildnbins, buildxmin, buildxmax);
    allhists.push_back(h);
    return h;
  } else {
    TFile *f = new TFile(histoutputfilename);
    auto h=(TH1F*)f->Get(hname);
    if (!h)
      cout << "histogram "<<hname<<" not found in file "<<histoutputfilename<<endl;
    else
      h->SetDirectory(0);
    f->Close();
    return h;
  }
}

TH1F *geth(TString hnametitle)
{
  return geth(hnametitle,hnametitle);
}

TH2F *geth2d(TString hname, TString htitle)
{
  if (firstRunMacro) {
    TH2F *h;
    h = new TH2F(hname+buildNamesuffix, buildTitlesuffix+htitle, buildnbins, buildxmin, buildxmax, buildnbinsy, buildymin, buildymax);
    allhists.push_back(h);
    return h;
  } else {
    TFile *f = new TFile(histoutputfilename);
    auto h=(TH2F*)f->Get(hname);
    if (!h)
      cout << "histogram "<<hname<<" not found in file "<<histoutputfilename<<endl;
    else
      h->SetDirectory(0);
    f->Close();
    return h;
  }
}

TH2F *geth2d(TString hnametitle)
{
  return geth2d(hnametitle,hnametitle);
}

void WriteAllHists()
{
  for (auto h:allhists)
    if (h!=0) h->Write();
}

void fitdphi(TH1F *h, float &sigma, float &error)
{
  //  mc_dphi
  //  TF1 *f = new TF1("myfit","exp((x-3.14216)/[0])/([0]*(1-exp(-3.1416/[0])))", 2./3*3.142, 3.1416);
  TF1 *f = new TF1("myfit","[1]*exp((x-3.14216)/[0])", 2./3*3.142, 3.1416);
  f->SetParameter(0,0.2);
  h->Fit(f,"NQ");

  cout<<"sigma = "<<f->GetParameter(0)<<endl;
  sigma = f->GetParameter(0);
  error = f->GetParError(0);
}



TString nicemeanstr(TH1F *h)
{
  float mean = h->GetMean();
  float std = h->GetMeanError();
  return TString::Format("%.3f #pm %.3f",round(mean*1000)/1000,round(std*1000)/1000);
}

TString nicewidthstr(TH1F *h)
{
  float s, e;
  fitdphi(h,s,e);
  return TString::Format("%.3f #pm %.3f",round(s*1000)/1000,round(e*1000)/1000);
}



TString plotsfolder = "plots";
bool PbPb = true;
TString plottitle = "";
TString plotfilenameend = "";
TString aktstring = "anti-k_{T}";
TString plotsecondline = "";
TString plotthirdline = "";
bool ploteffectiveentries = false;
bool plotcompatibility = false;

TString centralityLabel = "";


TString plotytitle = "Counts";
bool plotdivide = true;
bool plotputmean = false;
bool plotputwidth = false;
TString plotdatacaption = "data";
TString plotmccaption = "mc";
TString ploth11caption = "h11";
float plotymin = 9999;
float plotymin1 = 9999;
float plotymin2 = 9999;
float plotymax = 9999;
float plotymax1 = 9999;
float plotymax2 = 9999;
float plotyline = 9999;

//int ccounter = 0;
bool plotlegend = true;
LegendPosition plotlegendpos = TopRight;
bool plotylog = false;

float plottextposx = 0.55, plottextposy = 0.79;
float plotmeanposx = 0.2, plotmeanposy = 0.5;

bool plotoverwritecolors = true;

void Normalize(vector<TH1 *> hists)
{
  for (auto h:hists)
    h->Scale(1/h->Integral());
}

void NormalizeAllHists(vector<TH1 *> except=vector<TH1 *>(0))
{
  map<TString,int> namem;
  for (auto h:except) namem[h->GetName()] = 1;
  for (auto h:allhists)
    if (namem[h->GetName()]!=1)
      Normalize(allhists);
}

TLegend *getLegend()
{
  TLegend *l;
  if (plotlegendpos==TopRight)
    l = new TLegend(0.48,0.6,0.85,0.8);
  else if (plotlegendpos==BottomRight)
    l = new TLegend(0.48,0.2,0.85,0.4);
  else if (plotlegendpos==TopLeft)
    l = new TLegend(0.2,0.6,0.55,0.8);
  else if (plotlegendpos==BottomLeft)
    l = new TLegend(0.2,0.2,0.55,0.4);
  else 
    l= new TLegend(1,1,1,1);//outside

  //l->SetTextSize(20);

  return l;
}

TCanvas *Draw(vector<TH1 *> hists,TString options = "E1")
{
  ccounter++;
  TCanvas *c= new TCanvas(Form("%d",ccounter),Form("%d",ccounter),600,600);

  TLegend *l = getLegend();

  TLatex *Tl = new TLatex();



  int i=0;

  TString filename = "Draw_";

  for(auto h:hists) {
    filename+=h->GetName();
    if (plotoverwritecolors) {
      h->SetMarkerColor(TColor::GetColorDark(i+2));
      h->SetLineColor(TColor::GetColorDark(i+2));
    }
    l->AddEntry(h,h->GetTitle(),"P");
    if (i==0) {
      if (plotymax!=9999) h->SetMaximum(plotymax);
      if (plotymin!=9999) h->SetMinimum(plotymin);

      h->Draw(options);
    } else h->Draw(options+"same");
    i++;
    cout<<h->GetTitle()<<" : "<<h->GetEffectiveEntries()<<"("<<h->GetEntries()<<")"<<
      "mean : "<<h->GetMean()<<"±"<<h->GetMeanError()<<endl;
    if (plotputmean) {
      Tl->DrawLatexNDC(plotmeanposx, plotmeanposy-i*0.07, Form("%.3f#pm%.3f",h->GetMean(),h->GetMeanError()));
      TLine *l = new TLine(h->GetMean(),h->GetMinimum(),h->GetMean(),h->GetMaximum());
      l->SetLineColor(h->GetLineColor());
      l->SetLineWidth(2);
      l->Draw();

    }

  }
  cout<<endl;

 Tl->DrawLatexNDC(plottextposx, plottextposy, aktstring);
 if (plotsecondline!="")
   Tl->DrawLatexNDC(plottextposx, plottextposy-0.075,plotsecondline);
 if (plotthirdline!="")
   Tl->DrawLatexNDC(plottextposx, plottextposy-0.15,plotthirdline);

  if (plotyline!=9999) {
    TLine *line = new TLine(hists[0]->GetXaxis()->GetXmin(), plotyline, hists[0]->GetXaxis()->GetXmax(), plotyline);
    line->SetLineStyle(2);
    line->Draw();
  }



  if (plotlegend)
    l->Draw();
  if (plotylog)
    c->SetLogy();

  c->Update();

  gSystem->MakeDirectory("bjetplots");

  c->SaveAs(Form("bjetplots/%s.pdf",filename.Data()));
  return c;

}

void Print(TH1F *h)
{
    cout<<"Histogram "<<h->GetName()<<" : "<<h->GetTitle()<<endl;
  for (int i=1;i<=h->GetNbinsX();i++)
    cout<<"   "<<i<<" : "<<h->GetBinContent(i)<<" ± "<<h->GetBinError(i)<<endl;
}


void RenameBinLabelsX(TH1 *h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    h->GetXaxis()->SetBinLabel(i,Form("%d",i-1));
}

void RenameBinLabelsX(TH1 *h, vector<TString> labels)
{
  if (h->GetNbinsX()!=(int)labels.size()) cout<<" wrong number of labels/bins"<<endl;

  for (int i=1;i<=h->GetNbinsX();i++)
    h->GetXaxis()->SetBinLabel(i,labels[i-1]);
}

void RenameBinLabelsY(TH1 *h)
{
  for (int i=1;i<=h->GetNbinsY();i++)
    h->GetYaxis()->SetBinLabel(i,Form("%d",i-1));
}

void RenameBinLabelsY(TH1 *h, vector<TString> labels)
{
    if (h->GetNbinsY()!=(int)labels.size()) cout<<" wrong number of labels/bins"<<endl;
  for (int i=1;i<=h->GetNbinsY();i++)
    h->GetYaxis()->SetBinLabel(i,labels[i-1]);
}



map<TString, TString> histToStyle;

TString legendoption(TString drawoption)
{
  if (drawoption.Contains("hist")) return "L";
  else return "P";
}

void SetMC(vector<TH1F *> vh)
{
  for (auto h:vh) {
    histToStyle[h->GetName()] = "E1";//hist
    //h->SetLineWidth(2);
    //h->SetMarkerStyle(kNone);
    h->SetMarkerStyle(kOpenSquare);
  }
  //h->SetFillStyle(0);
}

void SetData(vector<TH1F *> vh)
{
  for (auto h:vh) {
    histToStyle[h->GetName()] = "E1";
    h->SetMarkerStyle(kFullCircle);
  }
}


void SetInc(vector<TH1F *> vh)
{
  for (auto h:vh) {
    h->SetLineColor(TColor::GetColorDark(kGray));
    h->SetMarkerColor(TColor::GetColorDark(kGray));
    h->SetFillColor(TColor::GetColorDark(kGray));
  }
}

void SetB(vector<TH1F *> vh)
{
  for (auto h:vh) {
    h->SetLineColor(darkred);
    h->SetMarkerColor(darkred);
    h->SetFillColor(darkred);
  }
}

void SetTruth(vector<TH1F *> vh)
{
  for (auto h:vh) {
    h->SetLineStyle(7);
    h->SetMarkerStyle(kOpenCircle);
  }
}


void DrawCompare(TH1F *h1, TH1F *h2, TString caption = "x_{J}",TH1F *h11=0)
{
  float textposx = 0.2, textposy = 0.77;
  TString title = plotytitle;
  float bw = h1->GetBinWidth(1);

  if (plotytitle=="Counts" && bw!=1)
    title = "Counts/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  TString legend3 = ploth11caption;
  //int color1 = kBlack;
  //int color2 = TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle()+h2->GetTitle(),title+h1->GetTitle()+h2->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (plotylog)
    pad1->SetLogy();

  float legendx1 = 0.55, legendy1 = 0.68, legendx2 = 0.84, legendy2 = 0.88;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetHeader(centralityLabel);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(gStyle->GetTextSize());

//legend1,2
 l->AddEntry(h1, h1->GetTitle(),legendoption(histToStyle[h1->GetName()]));
  if (h11!=0)
    l->AddEntry(h11, h11->GetTitle(),legendoption(histToStyle[h11->GetName()]));
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, h2->GetTitle(), legendoption(histToStyle[h2->GetName()]));

  if (plotymax!=9999) {
    h1->SetMaximum(plotymax);
    h2->SetMaximum(plotymax);
  }
  else {
    float ymax = plotylog ? pow(max(h1->GetMaximum(),h2->GetMaximum()),1.3):max(h1->GetMaximum(), h2->GetMaximum()) * 1.3;
    h1->SetMaximum(ymax);
    h2->SetMaximum(h1->GetMaximum());
  }
  if (plotymin!=9999) {
    //    h2->SetMaximum(h1->GetMaximum()*1.3);
    //    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(plotymin);
    h2->SetMinimum(plotymin);
  }


  h2->GetYaxis()->SetTitle(title);
  h2->GetXaxis()->SetLabelSize(0);

  //h1->SetMarkerColor(color1);
  //h1->SetLineColor(color1);
  //h2->SetMarkerColor(color2);
  //h2->SetLineColor(color2);
  //h2->SetFillColor(color2);
  //h2->SetFillStyle(1001);
  //if (h11!=0) {
  //h11->SetMarkerColor(darkblue);
  //h11->SetLineColor(darkblue);
  //}

  h2->Draw(histToStyle[h2->GetName()]);//"hist");
  h1->Draw(Form("%ssame",histToStyle[h1->GetName()].Data()));
  if (h11!=0) h11->Draw(Form("%ssame",histToStyle[h11->GetName()].Data()));
  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (plotsecondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.07,plotsecondline);
 if (plotthirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,plotthirdline);
 if (ploteffectiveentries)
   Tl->DrawLatexNDC(0.2,0.9,Form("Eff. entries data: %d, mc: %d",
         (int)h1->GetEffectiveEntries(),(int)h2->GetEffectiveEntries()));
 if (plottitle!="")
   Tl->DrawLatexNDC(0.2,0.9,plottitle);

  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  h3->SetMarkerColor(kGray+3);
  h3->SetLineColor(kGray+3);
  h3->SetLineStyle(1);
  h3->SetLineWidth(1);


  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(plotdivide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  //  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotdatacaption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotmccaption.Data())+nicemeanstr(h2));
    if (h11!=0)
      Tl2->DrawLatexNDC(textposx, 0.37, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),ploth11caption.Data())+nicemeanstr(h11));
    
  } else if (plotputwidth) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#sigma(%s)^{data} = ",caption.Data())+nicewidthstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#sigma(%s)^{mc} = ",caption.Data())+nicewidthstr(h2));

  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  TString h11title = h11==0?"":h11->GetTitle();
  TString outname = TString::Format("Compare%s_%s_%s%s%s",plottitle.Data(),h1->GetName(),h2->GetName(),h11title.Data(),plotfilenameend.Data());//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

  outname.ReplaceAll(" ","_");

  c1->SetName(outname);

  TString plotmacrofolder = plotsfolder+"/macro";
  gSystem->MakeDirectory(plotsfolder);
  gSystem->MakeDirectory(plotmacrofolder);  

  c1->SaveAs(Form("%s/%s.pdf",plotsfolder.Data(),outname.Data()));
  c1->SaveAs(Form("%s/%s.C",plotmacrofolder.Data(),outname.Data()));


}

void DrawCompare2(TH1F *h1, TH1F *h2, TString caption = "A_{J}",TH1F *h11=0)
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotytitle;
  float bw = h1->GetBinWidth(1);

  if (plotytitle=="Counts" && bw!=1)
    title = "Counts/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  TString legend3 = ploth11caption;
  int color1 = kBlack;
  int color2 = TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (plotylog)
    pad1->SetLogy();

  float legendx1 = 0.58, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetHeader(centralityLabel);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(h1->GetXaxis()->GetTitleSize());


  l->AddEntry(h1, legend1,"P");
  if (h11!=0)
    l->AddEntry(h11, legend3,"P");
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, legend2, "F");

  if (plotymax!=9999) {
    h1->SetMaximum(plotymax);
    h2->SetMaximum(plotymax);
  }
  else {
    float ymax = plotylog ? pow(max(h1->GetMaximum(),h2->GetMaximum()),1.3):max(h1->GetMaximum(), h2->GetMaximum()) * 1.3;
    h1->SetMaximum(ymax);
    h2->SetMaximum(h1->GetMaximum());
  }
  if (plotymin!=9999) {
    //    h2->SetMaximum(h1->GetMaximum()*1.3);
    //    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(plotymin);
    h2->SetMinimum(plotymin);
  }


  h2->GetYaxis()->SetTitle(title);
  h2->GetXaxis()->SetLabelSize(0);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);
  if (h11!=0) {
  h11->SetMarkerColor(darkblue);
  h11->SetLineColor(darkblue);
  }

  h2->Draw("hist");
  h1->Draw("same");
  if (h11!=0) h11->Draw("same");
  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (plotsecondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,plotsecondline);
 if (plotthirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,plotthirdline);
 if (ploteffectiveentries)
   Tl->DrawLatexNDC(0.2,0.9,Form("Eff. entries data: %d, mc: %d",
				 (int)h1->GetEffectiveEntries(),(int)h2->GetEffectiveEntries()));
 if (plottitle!="")
   Tl->DrawLatexNDC(0.2,0.9,plottitle);

  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(plotdivide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotdatacaption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotmccaption.Data())+nicemeanstr(h2));
    if (h11!=0)
      Tl2->DrawLatexNDC(textposx, 0.37, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),ploth11caption.Data())+nicemeanstr(h11));
    
  } else if (plotputwidth) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#sigma(%s)^{data} = ",caption.Data())+nicewidthstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#sigma(%s)^{mc} = ",caption.Data())+nicewidthstr(h2));

  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  TString h11title = h11==0?"":h11->GetTitle();

  c1->SaveAs(Form("%s/Compare%s_%s_%s%s.pdf",plotsfolder.Data(),plottitle.Data(),h1->GetTitle(),h2->GetTitle(),h11title.Data()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}

void DrawBoth(TH1F *h1, TH1F *h2, bool logy=false, TString caption = "A_{J}", TString h1title="", TString h2title="", TString secondline = "", TString thirdline = "")
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotytitle;
  float bw = h1->GetBinWidth(1);

  //  if (plotytitle=="Counts" && bw!=1)
  //    title = plotytitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";



  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  int color1 = kBlack;
  int color2 = kBlack;//TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (logy)
    pad1->SetLogy();

  float legendx1 = 0.58, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(h1->GetXaxis()->GetTitleSize());


  l->AddEntry(h1, legend1,"P");
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, legend2, "F");

  //  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  //  h2->SetMaximum(h1->GetMaximum());

  //  float ymax = plotylog ? pow(max(h1->GetMaximum(), h2->GetMaximum()), 1.3) : max(h1->GetMaximum(), h2->GetMaximum()) * 1.3;
  //  h1->SetMaximum(ymax);
  

  if (plotymin1!=9999) {
    h1->SetMinimum(plotymin1);
  }

  if (plotymin2!=9999) {
    h2->SetMinimum(plotymin2);
  }

  if (plotymax1!=9999) {
    h1->SetMaximum(plotymax1);
  }

  if (plotymin2!=9999) {
    h2->SetMaximum(plotymax2);
  }


  h2->GetYaxis()->SetTitle(h2title);
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetLabelSize(0);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);

  h2->Draw("E1");


  ///  h2->Draw("hist");
  ///  h1->Draw("same");
  //  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,thirdline);

  pad1->Draw();

  if (plotyline!=9999) {
    TLine *line = new TLine(gPad->GetUxmin(), plotyline, gPad->GetUxmax(), plotyline);
    line->SetLineStyle(2);
    line->Draw();
  }



  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.5);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  //  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  /*  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);
  */


  h1->GetYaxis()->SetTitle(h1title);//plotdivide ? "ratio" : "difference");
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitle(caption);
  h1->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h1->SetMarkerStyle(21);
  h1->Draw("E1");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotdatacaption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotmccaption.Data())+nicemeanstr(h2));
  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s.pdf",plotsfolder.Data(),h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}

vector<TString> stackCaptions = {};
vector<TH1F *> getHists(THStack *stack)
{
  vector<TH1F *> v;
  stackCaptions.clear();
  TList *lhist = stack->GetHists();
  if (lhist) {
    TH1F *h;
    TIter next(lhist);
    while ((h=(TH1F*)next())) {
      v.push_back(h);
      stackCaptions.push_back(h->GetTitle());
    }
  }
  return v;

}

TH1F *sumstack(THStack *s,TString name)
{
  auto hs = getHists(s);
  if (hs.size()==0) return 0;

  TH1F *h = (TH1F *)hs[0]->Clone(name);
  h->SetTitle(name);
  h->Reset();
  for (auto hi:hs)
    h->Add(hi);

  return h;
}

void normalizestack(THStack *s, float norm)
{
  auto hs = getHists(s);
  if (hs.size()==0) return;

  float sum = 0;
  for (auto hi:hs)
    sum+=hi->Integral();
  
   for (auto hi:hs) {
     hi->Scale(norm/sum);
   }
}




TCanvas * DrawStack(THStack *h, TH1F *hontop, TString xtitle, TString ytitle, float ymin = 999, float ymax = 999)
{
  auto c2 = getc();
  if (plotylog)
    c2->SetLogy();
  if (ymin!=999)
    h->SetMinimum(ymin);
  if (ymax!=999)
    h->SetMaximum(ymax);
  h->Draw("hist");
  if (hontop!=0)
    hontop->Draw("E1,same");

  TLegend *l = new TLegend(0.5, 0.7, 0.84, 0.84);
  auto hl = getHists(h);
  for (unsigned i=0;i<hl.size(); i++)
    l->AddEntry(hl[i], stackCaptions[i],"F");
  if (hontop!=0)
    l->AddEntry(hontop,"unmerged","P");
  l->Draw();

  //axes don't exist until Draw 
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  c2->Modified();
  c2->Update();
  c2->SaveAs(plotsfolder+"/"+h->GetTitle()+".pdf");

  return c2;
}

//awful...
TString getWithoutSuffix(TString s)
{
  if (s[s.Length()-1]==')' && s[s.Length()-3]=='(')
    return TString(s(0,s.Length()-3));

  return s;
}

THStack *stackhists(vector<TH1F *>hists, vector<int> color, TString name, TString suffix = "")
{
  THStack *hs = new THStack(name,name);
  int N = hists.size();
  for (int i=0;i<N;i++) {
  //for (int i=N-1;i>=0;i--) {
    hists[i]->SetFillColor(color[i]);
    hists[i]->SetLineColor(color[i]);
    hists[i]->SetMarkerColor(color[i]);
    hists[i]->SetFillStyle(1001);
    hists[i]->SetTitle(getWithoutSuffix(hists[i]->GetTitle())+suffix);

    hs->Add(hists[i],"hist");
  }

  return hs;
}


void DrawCompare(TH1F *h1, THStack *hstack, TString caption = "x_{J}",TString stupidname = "")
{
  float textposx = 0.2, textposy = 0.77;
  TString title = plotytitle;

  //auto h2 = sumstack(hstack, Form("%ssum",hstack->GetTitle()));
  normalizestack(hstack, h1->Integral());

  auto h2 = sumstack(hstack, Form("%ssum%s",hstack->GetName(),stupidname.Data()));
  h2->ResetStats();

  float bw = h2->GetBinWidth(2);

  if (plotytitle=="Counts" && bw!=1)
    title = plotytitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  int color1 = kBlack;
  int color2 = TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (plotylog)
    pad1->SetLogy();

  float legendx1 = 0.55, legendy1 = 0.68, legendx2 = 0.84, legendy2 = 0.88;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetHeader(centralityLabel);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(gStyle->GetTextSize());

  l->AddEntry(h1, legend1,"P");

  auto hhh = getHists(hstack);
  // if (hhh.size()==3) {
  //   l->AddEntry(hhh[0], "B", "F");
  //   l->AddEntry(hhh[1], "C", "F");
  //   l->AddEntry(hhh[2], "L", "F");
  // }else
  // if (hhh.size()==4) {
  //   l->SetY1(0.58);
  //   l->AddEntry(hhh[0], "B primary", "F");
  //   l->AddEntry(hhh[1], "B GSP", "F");
  //   l->AddEntry(hhh[2], "C", "F");
  //   l->AddEntry(hhh[3], "L", "F");
  // }else
  //   l->AddEntry(hstack, "", "F");
  
  for (int i=hhh.size()-1;i>=0;i--)
    l->AddEntry(hhh[i],hhh[i]->GetTitle(),"F");


  // float ymax = plotylog ? pow(max(h1->GetMaximum(), h2->GetMaximum()), 1.3) : max(h1->GetMaximum(), h2->GetMaximum()) * 1.3;
  // h1->SetMaximum(ymax);
  // h2->SetMaximum(h1->GetMaximum());

  // //  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  // //  h2->SetMaximum(h1->GetMaximum());

  // if (plotymin!=9999) {
  //   h2->SetMaximum(h1->GetMaximum()*1.3);
  //   h1->SetMaximum(h1->GetMaximum()*1.3);

  //   hstack->SetMinimum(plotymin);
  //   h1->SetMinimum(plotymin);
  //   h2->SetMinimum(plotymin);
  // }


 if (plotymax!=9999) {
    h1->SetMaximum(plotymax);
    h2->SetMaximum(plotymax);
    hstack->SetMaximum(plotymax);
  }
  else {
    float ymax = plotylog ? pow(max(h1->GetMaximum(),h2->GetMaximum()),1.3):max(h1->GetMaximum(), h2->GetMaximum()) * 1.3;
    h1->SetMaximum(ymax);
    h2->SetMaximum(ymax);
    hstack->SetMaximum(ymax);
  }
  if (plotymin!=9999) {
    //    h2->SetMaximum(h1->GetMaximum()*1.3);
    //    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(plotymin);
    h2->SetMinimum(plotymin);
        hstack->SetMinimum(plotymin);
  }


  h2->GetYaxis()->SetTitle(title);
  h2->GetXaxis()->SetLabelSize(0);


  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  //  h2->SetMarkerColor(color2);
  //  h2->SetLineColor(color2);
  //  h2->SetFillColor(color2);
  //  h2->SetFillStyle(1001);

  hstack->Draw("hist");
  h1->Draw("P,same");


  hstack->GetXaxis()->SetLabelSize(0);

  //h2->Draw();
  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (plotsecondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,plotsecondline);
 if (plotthirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,plotthirdline);
 if (ploteffectiveentries)
   Tl->DrawLatexNDC(0.2,0.9,Form("Eff. entries data: %d, mc: %d",
				 (int)h1->GetEffectiveEntries(),(int)h2->GetEffectiveEntries()));
 if (plotcompatibility)
  Tl->DrawLatexNDC(0.2,0.9,Form("#chi^2 : %.5f, KS : %.5f ",
    h1->Chi2Test(h2, "UW"),h1->KolmogorovTest(h2)));  

  cout<<"Compatibility : chi2 = "<<h1->Chi2Test(h2, "UW")<<", KS = "<<h1->KolmogorovTest(h2)<<endl;

  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(plotdivide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotdatacaption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{%s} = ",caption.Data(),plotmccaption.Data())+nicemeanstr(h2));
  } else if (plotputwidth) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#sigma(%s)^{%s} = ",caption.Data(),plotdatacaption.Data())+nicewidthstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#sigma(%s)^{%s} = ",caption.Data(),plotmccaption.Data())+nicewidthstr(h2));

  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s_%s.pdf",plotsfolder.Data(),h1->GetTitle(),h2->GetTitle()));

}




class macro
{
  TString macroname;
  bool firstRun;
  TString filename;
public:
  TString getname() { return macroname;}
  macro(TString name, bool firstrun = true) 
  {
    firstRun = firstrun;
    macroname = name;
    filename = Form("%s_hists.root",macroname.Data());
    firstRunMacro = firstrun;
    histoutputfilename = filename;
    cout<<"Starting macro "<<macroname<<"!"<<endl;
  }
  ~macro() 
  {
    if (!firstRun) return;
    cout<<"Writing histograms to the file "<<filename<<endl;
    TFile *f = new TFile(filename,"recreate");
    for (auto x:allhists)
      x->Write();
    f->Close();
    cout<<"done!"<<endl;
  }
};

// macro StartMacro(TString macroName)
// {
//   macro m(macroName);
//   return m;
// }


#endif

