#include "./CMS_lumi.C"
#include "../helpers/plotting.h"

void plotDphi(int inc_or_bjet=0)
{
  macro m("plotDphi");

 writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  lumi_sqrtS = "404 #mub^{-1} (5.02 TeV PbPb)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  int iPos=11; // : top-left, left-aligned
  //int  iPos=33;// : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  
  TFile *fin = new TFile("results_0702_newweight/xJdphi.root");

  string species = "inc";
  if(inc_or_bjet) species = "bjt";

  TH1F *hData010 = (TH1F*) fin->Get(Form("dphi_data_%s_0_10",species.c_str()));
  TH1F *hData1030 = (TH1F*) fin->Get(Form("dphi_data_%s_10_30",species.c_str()));
  TH1F *hData30100 = (TH1F*) fin->Get(Form("dphi_data_%s_30_100",species.c_str()));
  TH1F *hDataPP = (TH1F*) fin->Get(Form("dphi_data_%s_pp",species.c_str()));
  
  TH1F *hMC010 = (TH1F*) fin->Get(Form("dphi_mc_%s_0_10",species.c_str()));
  TH1F *hMC1030 = (TH1F*) fin->Get(Form("dphi_mc_%s_10_30",species.c_str()));
  TH1F *hMC30100 = (TH1F*) fin->Get(Form("dphi_mc_%s_30_100",species.c_str()));
  TH1F *hMCPP = (TH1F*) fin->Get(Form("dphi_mc_%s_pp",species.c_str()));

  TF1 *fData010 = (TF1*) fin->Get(Form("fit_dphi_data_%s_0_10",species.c_str()));
  TF1 *fData1030 = (TF1*) fin->Get(Form("fit_dphi_data_%s_10_30",species.c_str()));
  TF1 *fData30100 = (TF1*) fin->Get(Form("fit_dphi_data_%s_30_100",species.c_str()));
  TF1 *fDataPP = (TF1*) fin->Get(Form("fit_dphi_data_%s_pp",species.c_str()));
  
  TF1 *fMC010 = (TF1*) fin->Get(Form("fit_dphi_mc_%s_0_10",species.c_str()));
  TF1 *fMC1030 = (TF1*) fin->Get(Form("fit_dphi_mc_%s_10_30",species.c_str()));
  TF1 *fMC30100 = (TF1*) fin->Get(Form("fit_dphi_mc_%s_30_100",species.c_str()));
  TF1 *fMCPP = (TF1*) fin->Get(Form("fit_dphi_mc_%s_pp",species.c_str()));



  // if(inc_or_bjet) hMCPP = (TH1F*) fin->Get(Form("dphi_data_%s_0_10;1",species.c_str()));

  TCanvas *c1=new TCanvas("c1","c1",600,600);

  int color = kblue;
  if(inc_or_bjet) color=kred;

  
  hData010->SetMarkerColor(color);
  hData010->SetLineColor(color);
  hMC010->SetLineColor(color);

  hData010->GetXaxis()->CenterTitle(1);
  hData010->GetYaxis()->CenterTitle(1);
  hData010->SetYTitle("Event fraction");
  hData010->SetXTitle("|#Delta#phi|");
  hData010->SetMinimum(0.001);
  hData010->Draw();
  hMC010->SetMinimum(0.001);
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

 SavePlot(c1,"dphi010"+species);

  TCanvas *c2=new TCanvas("c2","c2",600,600);


  hData1030->SetMarkerColor(color);
  hData1030->SetLineColor(color);
  hMC1030->SetLineColor(color);

  hData1030->GetXaxis()->CenterTitle(1);
  hData1030->GetYaxis()->CenterTitle(1);
  hData1030->SetYTitle("Event fraction");
  hData1030->SetXTitle("|#Delta#phi|");
  hData1030->Draw();
  hMC1030->SetMarkerSize(0);
  hMC1030->Draw("h,same");
  
  
  CMS_lumi(c2, iPeriod, iPos );

  l->Draw();

 SavePlot(c2,"dphi1030"+species);

    TCanvas *c3=new TCanvas("c3","c3",600,600);


  hData30100->SetMarkerColor(color);
  hData30100->SetLineColor(color);
  hMC30100->SetLineColor(color);

  hData30100->GetXaxis()->CenterTitle(1);
  hData30100->GetYaxis()->CenterTitle(1);
  hData30100->SetYTitle("Event fraction");
  hData30100->SetXTitle("|#Delta#phi|");
  hData30100->Draw();
  hMC30100->SetMarkerSize(0);
  hMC30100->Draw("h,same");
  
  
  CMS_lumi(c3, iPeriod, iPos );

  lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp)";

    l->Draw();

SavePlot(c3,"dphi30100"+species);


   TCanvas *c4=new TCanvas("c4","c4",600,600);


  hDataPP->SetMarkerColor(color);
  hDataPP->SetLineColor(color);
  hMCPP->SetLineColor(color);

  hDataPP->GetXaxis()->CenterTitle(1);
  hDataPP->GetYaxis()->CenterTitle(1);
  hDataPP->SetYTitle("Event fraction");
  hDataPP->SetXTitle("|#Delta#phi|");
  hDataPP->Draw();
  hMCPP->SetMarkerSize(0);
  hMCPP->Draw("h,same");
  
  
  CMS_lumi(c4, iPeriod, iPos );


  
  l->Draw();
SavePlot(c4,"dphipp"+species);

 

cout<< " & Pythia & Data  \\\\ "<<endl;
cout<<"\\hline"<<endl;
cout << std::fixed <<setprecision(3);
cout<<"pp & "        <<fMCPP->GetParameter(0)   <<"\\pm"<<fMCPP->GetParError(0) <<" & "<<fDataPP->GetParameter(0)<<"\\pm"<<fDataPP->GetParError(0)<<" \\\\"<<endl;
cout<<"30-100\\% & " <<fMC30100->GetParameter(0)<<"\\pm"<<fMC30100->GetParError(0) <<" & "<<fData30100->GetParameter(0)<<"\\pm"<<fData30100->GetParError(0)<<" \\\\"<<endl;
cout<<"10-30\\% & "  <<fMC1030->GetParameter(0) <<"\\pm"<<fMC1030->GetParError(0) <<" & "<<fData1030->GetParameter(0)<<"\\pm"<<fData1030->GetParError(0)<<" \\\\"<<endl;
cout<<"0-10\\% & "   <<fMC010->GetParameter(0)  <<"\\pm"<<fMC010->GetParError(0) <<" & "<<fData010->GetParameter(0)<<"\\pm"<<fData010->GetParError(0)<<" \\\\"<<endl;



}
