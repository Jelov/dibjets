#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"

vector<float> bins = {0,20,60,200};
vector<TString> binnames = {"0-10%", "10-30%", "30-100%"}; 

int getbinindex(float bin)
{
  for(unsigned i=0;i<bins.size();i++)
    if (bins[i]>bin) return i-1;
  return bins.size()-2;
}


float weight1SL(dict d)
{
  float w = d["weight"];
  if (d["pairCodeSL1"]==0) w*=processweight((int)d["bProdCode"]);
  return w;
}

bool IsSignal(dict d) { return d["subidSL"]==0 && d["refptSL"]>20;}

bool NearSide(dict d)
{
  return d["dphiSL1"]<PI23;
}

bool AwaySide(dict d)
{
  return d["dphiSL1"]>PI23;
}

void dphisignalconsistency()
{

  auto fmcpp = config.getfile_djt("mcppbfa");
  auto fmcPb = config.getfile_djt("mcPbbfa");

  buildh(3,0,PI);
  auto hmcppdphi = geth("hmcppdphi","MC pp;Signal only #Delta#phi");
  buildh(1,0,PI13);
  auto hmcppdphiNS = geth("hmcppdphiNS","MC pp NS;#Delta#phi");
  buildh(1,PI23,PI);
  auto hmcppdphiAS = geth("hmcppdphiAS","MC pp AS;#Delta#phi");


  int Nb = bins.size()-1;

  vector<TH1 *>hmcPbdphi(Nb),hmcPbdphiNS(Nb),hmcPbdphiAS(Nb),hmcPbdphiNSall(Nb),hmcPbdphiASall(Nb),hmcPbdphibkg(Nb);
  for (int i=0;i<Nb;i++) {
    buildh(3,0,PI);
    hmcPbdphi[i] = geth(Form("hmcPbdphi%d",i),Form("MC Pb %s;Signal only #Delta#phi",binnames[i].Data()));
    hmcPbdphibkg[i] = geth(Form("hmcPbdphibkg%d",i),Form("MC Pb %s;Hydjet only #Delta#phi",binnames[i].Data()));
    buildh(1,0,PI13);
    hmcPbdphiNS[i] = geth(Form("hmcPbdphiNS%d",i),Form("MC NS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSall[i] = geth(Form("hmcPbdphiNSall%d",i),Form("MC NS ALL Pb %s;#Delta#phi",binnames[i].Data()));
    buildh(1,PI23,PI);
    hmcPbdphiAS[i] = geth(Form("hmcPbdphiAS%d",i),Form("MC AS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiASall[i] = geth(Form("hmcPbdphiASall%d",i),Form("MC AS ALL Pb %s;#Delta#phi",binnames[i].Data()));
  }

  Fill(fmcpp,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","dphiSL1","refparton_flavorForB1","pairCodeSL1","discr_csvV1_1"},[&] (dict d) {
    float w = d["weight"];
    if (d["pthat"]<80) return;
    if (d["refpt1"]<50)  return;

    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) {//&& d["discr_csvV1_1"]>0.9
      hmcppdphi->Fill(d["dphiSL1"],weight1SL(d));
      if (NearSide(d)) hmcppdphiNS->Fill(d["dphiSL1"],weight1SL(d));
      if (AwaySide(d)) hmcppdphiAS->Fill(d["dphiSL1"],weight1SL(d));
    }

  });//process inly 20% of data

  Fill(fmcPb,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","refptSL","dphiSL1","refparton_flavorForB1","subidSL","bin","pairCodeSL1","discr_csvV1_1"},[&] (dict d) {
    float w = d["weight"];
    if (d["pthat"]<80) return;
    if (d["refpt1"]<50) return;

    int b = getbinindex(d["bin"]);
    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut && IsSignal(d)) { //d["discr_csvV1_1"]>0.9
      hmcPbdphi[b]->Fill(d["dphiSL1"],weight1SL(d));
      if (NearSide(d)) hmcPbdphiNS[b]->Fill(d["dphiSL1"],weight1SL(d));
      if (AwaySide(d)) hmcPbdphiAS[b]->Fill(d["dphiSL1"],weight1SL(d));
    }
    
    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) { //d["discr_csvV1_1"]>0.9
      if (NearSide(d)) hmcPbdphiNSall[b]->Fill(d["dphiSL1"],weight1SL(d));
      if (AwaySide(d)) hmcPbdphiASall[b]->Fill(d["dphiSL1"],weight1SL(d));
    }

    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut && !IsSignal(d))
      hmcPbdphibkg[b]->Fill(d["dphiSL1"],weight1SL(d));

  });

  vector<float> x, y;
  for (int i=0;i<Nb;i++) {
    double pp1 = hmcppdphiNS->GetBinContent(1)*1E9;
    double pp3 = hmcppdphiAS->GetBinContent(1)*1E9;
    double PbPb1 = hmcPbdphiNSall[i]->GetBinContent(1)*1E9;
    double PbPb3 = hmcPbdphiASall[i]->GetBinContent(1)*1E9;

    float y_ = (PbPb3-PbPb1)/(pp3-pp1);
    float x_ = 0.5*(PbPb3+PbPb1-y_*(pp3+pp1))*1E-9;
    y.push_back(y_);
    x.push_back(x_);

    cout<<"bin "<<binnames[i]<<" x*1E9 = "<<x_*1E9<<" ; y = "<<y_<<endl;
    cout<<" bkg fraction in NS = "<<x_*1E9/PbPb1<<endl;
  }


  plotlegendpos = TopLeft;
  plotymin = 0; plotymax = 1E-7;
  plotyline = x[0];
  Draw({hmcPbdphibkg[0]});
  //(new TLine(hmcPbdphibkg[0]->GetXaxis()->GetXmin(),x[0],hmcPbdphibkg[0]->GetXaxis()->GetXmax(),x[0]))->Draw();

  plotymin = 0; plotymax = 6E-9;
  plotyline = x[1];
  Draw({hmcPbdphibkg[1]});
  //(new TLine(hmcPbdphibkg[1]->GetXaxis()->GetXmin(),x[1],hmcPbdphibkg[1]->GetXaxis()->GetXmax(),x[1]))->Draw();

  plotymin = -0.5E-9; plotymax = 2E-9;
  plotyline = x[2];
  Draw({hmcPbdphibkg[2]});
  //(new TLine(hmcPbdphibkg[2]->GetXaxis()->GetXmin(),x[2],hmcPbdphibkg[2]->GetXaxis()->GetXmax(),x[2]))->Draw();

  plotyline = 9999;
  plotymax = 1;

  hmcppdphiNS->Divide(hmcppdphiAS);
  for (int i=0;i<Nb;i++)
    hmcPbdphiNS[i]->Divide(hmcPbdphiAS[i]);

  buildh(Nb+1,0,Nb+1);
  auto hNSASratio = geth("hNSASratio","ratio of NS to AS");
  hNSASratio->SetBinContent(1,hmcppdphiNS->GetBinContent(1));
  hNSASratio->SetBinError(1,hmcppdphiNS->GetBinError(1));
  for (int i=0;i<Nb;i++) {
    hNSASratio->SetBinContent(i+2,hmcPbdphiNS[i]->GetBinContent(1));
    hNSASratio->SetBinError(i+2,hmcPbdphiNS[i]->GetBinError(1));
  }


  vector<TString> axnames;
  axnames.push_back("pp");
  for (auto s:binnames) axnames.push_back(s);
    RenameBinLabelsX(hNSASratio,axnames);


  hmcPbdphi.push_back(hmcppdphi);
  Normalize(hmcPbdphi);

  aktstring = "";
  plotlegendpos = TopLeft;
  Draw(hmcPbdphi);

  plotymin = 0;
  plotymax = 0.15;
  Draw({hNSASratio});

}

void hydjetestimation()
{
  dphisignalconsistency();
}