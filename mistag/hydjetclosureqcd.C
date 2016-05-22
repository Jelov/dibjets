#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

bool IsSignal2(dict d)
{
  return d["subid2"]==0 && d["refpt2"]>20;
}

bool NearSide2(dict d)
{
  //return d["dphiSL1"]<PI13;

  float dphi = d["dphi21"];
  float deta = abs(d["jteta1"]-d["jteta2"]); 
  return (dphi<PI13 && (dphi*dphi+deta*deta)>1) || (dphi>PI13 && ((dphi-PI13)*(dphi-PI13)+deta*deta)<1);
}

bool AwaySide2(dict d)
{
  return d["dphi21"]>PI23;
}

void checkclosure()
{
  vector<TH1F *>hsig(Nbins);
  vector<TH1F *>hasd(Nbins);
  vector<TH1F *>hbkg(Nbins);
  vector<TH1F *>hsub(Nbins);
  vector<TH1F *>hhyj(Nbins);
  vector<TH1F *>hshj(Nbins);
  vector<TH1F *>hsbn(Nbins);


  for (int i=0;i<Nbins;i++) {
    buildh(10,0,1);
    hsig[i] = geth(Form("hsig%d",i),Form("Signal away-side %s;x_{J}",binnames[i].Data())) ;
    hasd[i] = geth(Form("hasd%d",i),Form("Measured away-side %s;x_{J}",binnames[i].Data()));
    hbkg[i] = geth(Form("hbkg%d",i),Form("Near-side %s;x_{J}",binnames[i].Data()));
    hhyj[i] = geth(Form("hhyj%d",i),Form("Near-side hydjet %s;x_{J}",binnames[i].Data()));
    hsub[i] = geth(Form("hsub%d",i),Form("Subtracted NS %s;x_{J}",binnames[i].Data()));
    hshj[i] = geth(Form("hshj%d",i),Form("Subtracted Hydjet %s;x_{J}",binnames[i].Data()));
    hsbn[i] = geth(Form("hsbn%d",i),Form("Subtracted Naive %s;x_{J}",binnames[i].Data()));
  }




  auto fmcPb = config.getfile_djt("mcPbqcd");

  Fill(fmcPb,{"pthat","weight","jtpt1","refpt1","jtpt2","dphi21","subid2","refpt2","jteta1","jteta2","bin"},[&] (dict d) {
      if (d["pthat"]<pthatcut) return;
      
      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && d["jtpt2"]>pt2cut) {
        int bin = getbinindex(d["bin"]);
        
        float xj = d["jtpt2"]/d["jtpt1"];
        float w = d["weight"];
        if (AwaySide2(d)) hasd[bin]->Fill(xj, w);
        if (AwaySide2(d) && IsSignal2(d)) hsig[bin]->Fill(xj,w);

        if (NearSide2(d)) hbkg[bin]->Fill(xj,w);
        if (NearSide2(d) && !IsSignal2(d)) hhyj[bin]->Fill(xj,w);
      }
        



      });



  for (int i=0;i<Nbins;i++) {
    hsub[i]->Add(hasd[i],hbkg[i],1,-1);
    hsbn[i]->Add(hasd[i],hbkg[i],1,-1);
    hshj[i]->Add(hasd[i],hhyj[i],1,-1);
  }
//  for (int i=0;i<Nbins;i++) 
//    hincsub[i]->Add(hincasd[i],hincbkg[i],1,-1);

  buildh(bins);//Nbins,0,100);
  auto hcentrSubSIG = geth("hcentrSubSIG","Signal;bin;<x_{J}>");
  auto hcentrSubASD = geth("hcentrSubASD","Unsubtracted;bin;<x_{J}>");

  auto hcentrSubBKS = geth("hcentrSubBKS","Naive subtraction;bin;<x_{J}>");
  auto hcentrSubCLS = geth("hcentrSubCLS","Subtracted;bin;<x_{J}>");
  auto hcentrSubHJS = geth("hcentrSubHJS","Subtracted Hydjet;bin;<x_{J}>");



  for (int i=0;i<Nbins;i++) {
    hcentrSubSIG->SetBinContent(i+1,hsig[i]->GetMean());hcentrSubSIG->SetBinError(i+1,hsig[i]->GetMeanError());
    hcentrSubASD->SetBinContent(i+1,hasd[i]->GetMean());hcentrSubASD->SetBinError(i+1,hasd[i]->GetMeanError());
    hcentrSubBKS->SetBinContent(i+1,hsbn[i]->GetMean());hcentrSubBKS->SetBinError(i+1,hsbn[i]->GetMeanError());

    hcentrSubCLS->SetBinContent(i+1,hsub[i]->GetMean());hcentrSubCLS->SetBinError(i+1,hsub[i]->GetMeanError());
    hcentrSubHJS->SetBinContent(i+1,hshj[i]->GetMean());hcentrSubHJS->SetBinError(i+1,hshj[i]->GetMeanError());

    Draw({hsig[i],hsub[i],hshj[i]});
  }


  plotymin = 0.65;//0.4;
  plotymax = 0.75;//0.8;
  plotlegendpos = BottomRight;
  aktstring = "";


  plotputmean = false;
  //hcentrSubHJS - hydjet only subtraction
  SetMC({hcentrSubSIG, hcentrSubBKS, hcentrSubASD});
  SetData({hcentrSubCLS});
  hcentrSubSIG->SetMarkerColor(TColor::GetColorDark(2)); hcentrSubSIG->SetLineColor(TColor::GetColorDark(2));
  hcentrSubBKS->SetMarkerColor(TColor::GetColorDark(3)); hcentrSubBKS->SetLineColor(TColor::GetColorDark(3));
  hcentrSubASD->SetMarkerColor(TColor::GetColorDark(4)); hcentrSubASD->SetLineColor(TColor::GetColorDark(4));
  hcentrSubCLS->SetMarkerColor(TColor::GetColorDark(3)); hcentrSubCLS->SetLineColor(TColor::GetColorDark(3));

  plotoverwritecolors = false;
  Draw({hcentrSubSIG, hcentrSubBKS, hcentrSubASD, hcentrSubCLS});


  auto syst = (TH1F *)hcentrSubSIG->Clone("syst");
  syst->Add(hcentrSubCLS,-1);
  map<TString,float> m;
  for (unsigned i=0;i<bins.size()-1;i++)
    m[Form("closure%d%d",(int)bins[i],(int)bins[i+1])]=syst->GetBinContent(i+1);

  WriteToFile(plotfoldername+"/hydjetclosureqcdsyst.root",m);

}



void hydjetclosureqcd(bool firstRun = true)
{
  macro m("hydjetclosureqcd",firstRun);
  //looptupledryrun = true;
  checkclosure();
}