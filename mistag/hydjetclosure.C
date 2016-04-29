#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

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




  auto fmcPb = config.getfile_djt("mcPbbfa");

  Fill(fmcPb,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","refptSL","dphiSL1","refparton_flavorForB1","subidSL","bin","pairCodeSL1","discr_csvV1_1","jteta1","jtetaSL"},[&] (dict d) {
      if (d["pthat"]<80) return;
      
      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) {
        int bin = getbinindex(d["bin"]);
        
        float xj = d["jtptSL"]/d["jtpt1"];
        float w = weight1SL(d);
        if (AwaySide(d)) hasd[bin]->Fill(xj, w);
        if (AwaySide(d) && IsSignal(d)) hsig[bin]->Fill(xj,w);

        if (NearSide(d)) hbkg[bin]->Fill(xj,w);
        if (NearSide(d) && !IsSignal(d)) hhyj[bin]->Fill(xj,w);
      }
        



      });



  for (int i=0;i<Nbins;i++) {
    hsub[i]->Add(hasd[i],hbkg[i],1,-1*bkgfractionInNearSide[i]);
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
  }


  plotymin = 0.6;//0.4;
  plotymax = 0.75;//0.8;
  plotlegendpos = BottomRight;
  aktstring = "";


  plotputmean = false;
  //hcentrSubHJS - hydjet only subtraction
  SetMC({hcentrSubSIG, hcentrSubBKS, hcentrSubASD, hcentrSubCLS});
  Draw({hcentrSubSIG, hcentrSubBKS, hcentrSubASD, hcentrSubCLS});



}



void hydjetclosure()
{
  macro m("hydjetclosure");
  //looptupledryrun = true;
  checkclosure();
}