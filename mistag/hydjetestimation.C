#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

void dphisignalconsistency()
{

  auto fmcpp = config.getfile_djt("mcppbfa");
  auto fmcPb = config.getfile_djt("mcPbbfa");

  buildh(3,0,PI);
  auto hmcppdphi = geth("hmcppdphi","MC pp;Signal only #Delta#phi");
  buildh(1,0,PI);
  auto hmcppdphiNS = geth("hmcppdphiNS","MC pp NS;#Delta#phi");
  auto hmcppdphiAS = geth("hmcppdphiAS","MC pp AS;#Delta#phi");


  int Nb = bins.size()-1;

  vector<TH1 *> hmcPbdphi(Nb),hmcPbdphiNS(Nb),hmcPbdphiAS(Nb),hmcPbdphiNSall(Nb),hmcPbdphiASall(Nb),
                hmcPbdphibkg(Nb),hmcPbdphiNSbkg(Nb),hmcPbdphiASbkg(Nb),hmcPbdphiNSASratio(Nb);
                
  for (int i=0;i<Nb;i++) {
    buildh(3,0,PI);
    hmcPbdphi[i] = geth(Form("hmcPbdphi%d",i),Form("MC Pb %s;Signal only #Delta#phi",binnames[i].Data()));
    hmcPbdphibkg[i] = geth(Form("hmcPbdphibkg%d",i),Form("MC Pb %s;Hydjet only #Delta#phi",binnames[i].Data()));
    buildh(1,0,PI);
    hmcPbdphiNS[i] = geth(Form("hmcPbdphiNS%d",i),Form("MC NS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSbkg[i] = geth(Form("hmcPbdphiNSbkg%d",i),Form("MC NS Background Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSall[i] = geth(Form("hmcPbdphiNSall%d",i),Form("MC NS ALL Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiAS[i] = geth(Form("hmcPbdphiAS%d",i),Form("MC AS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiASbkg[i] = geth(Form("hmcPbdphiASbkg%d",i),Form("MC AS Background Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiASall[i] = geth(Form("hmcPbdphiASall%d",i),Form("MC AS ALL Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSASratio[i] = geth(Form("hmcPbdphiNSASratio%d",i),Form("MC NS AS ratio Pb %s;#Delta#phi",binnames[i].Data()));
  }

  Fill(fmcpp,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","dphiSL1","refparton_flavorForB1","pairCodeSL1","discr_csvV1_1","jteta1","jtetaSL"},[&] (dict d) {
    float w = weight1SLpp(d);
    if (d["pthat"]<pthatcut) return;
    if (d["refpt1"]<50)  return;

    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) {//&& d["discr_csvV1_1"]>0.9
      hmcppdphi->Fill(d["dphiSL1"],w);
      if (NearSide(d)) hmcppdphiNS->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hmcppdphiAS->Fill(d["dphiSL1"],w);
    }

  });

  Fill(fmcPb,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","refptSL","dphiSL1","refparton_flavorForB1","subidSL","bin","pairCodeSL1","discr_csvV1_1","jteta1","jtetaSL"},[&] (dict d) {
    float w = weight1SLPbPb(d);
    if (d["pthat"]<pthatcut) return;
    if (d["refpt1"]<50) return;

    int b = getbinindex(d["bin"]);
    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut && IsSignal(d)) { //d["discr_csvV1_1"]>0.9
      hmcPbdphi[b]->Fill(d["dphiSL1"],w);
      if (NearSide(d)) hmcPbdphiNS[b]->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hmcPbdphiAS[b]->Fill(d["dphiSL1"],w);
    }
    
    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) { //d["discr_csvV1_1"]>0.9
      if (NearSide(d)) hmcPbdphiNSall[b]->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hmcPbdphiASall[b]->Fill(d["dphiSL1"],w);
    }

    if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut && !IsSignal(d)) {
      hmcPbdphibkg[b]->Fill(d["dphiSL1"],w);
      if (NearSide(d)) hmcPbdphiNSbkg[b]->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hmcPbdphiASbkg[b]->Fill(d["dphiSL1"],w);
    }

  });

  vector<float> x;//, y;
  vector<TH1F *>vhx;
  for (int i=0;i<Nb;i++) {
    float pp1 = hmcppdphiNS->GetBinContent(1);
    float pp3 = hmcppdphiAS->GetBinContent(1);
    //cross-check: gives correct 0.755140, 0.261050, 0.0134587
    // float pp1 = hmcPbdphiNS[i]->GetBinContent(1);
    // float pp3 = hmcPbdphiAS[i]->GetBinContent(1);
    
    float PbPb1 = hmcPbdphiNSall[i]->GetBinContent(1);
    float PbPb3 = hmcPbdphiASall[i]->GetBinContent(1);

    cout<<pp1*1E9<<" "<<pp3*1E9<<" "<<PbPb1*1E9<<" "<<PbPb3*1E9<<endl;


    buildh(1,0,PI);
    auto hB = geth(Form("hB%d",i));
    auto hxx = geth(Form("hxx%d",i),Form("Estimated background in %s",binnames[i].Data()));
    hB->Divide(hmcppdphiNS,hmcppdphiAS);
    hxx->Divide(hmcPbdphiASall[i],hmcPbdphiNSall[i]);
    hxx->Multiply(hB);
    hxx->SetBinContent(1,1-hxx->GetBinContent(1));
    hB->SetBinContent(1,1-hB->GetBinContent(1));
    hxx->Divide(hB);


    float B = pp1/pp3;
    float xx_ = (1-PbPb3/PbPb1*B)/(1-B);

    cout<<" checking hist arithmetics "<<hxx->GetBinContent(1)<<" ± "<<hxx->GetBinError(1)<<" =? "<<xx_<<endl;


    hxx->Multiply(hmcPbdphiNSall[i]);
    vhx.push_back(hxx);


    //just for fun... turn on quenching of the away-side by alpha
    double alpha = 0.5;
    double Bprime = 1/alpha*B;
    float xalpha_ = (1-PbPb3/PbPb1*Bprime)/(1-Bprime);
    cout<<" bkg fraction in NS = "<<xx_<<", quenched by "<<alpha<<" bkg = "<<xalpha_<<endl;
    cout<<setprecision(10)<<xx_*PbPb1<<endl;

    x.push_back(xx_*PbPb1);

    bkgfractionInNearSide[i] = xx_;
  }



  plotlegendpos = TopLeft;
  plotymin = 0; plotymax = 1E-7;
  plotyline = x[0];
  Draw({hmcPbdphibkg[0],vhx[0]});

  plotymin = 0; plotymax = 1E-8;
  plotyline = x[1];
  Draw({hmcPbdphibkg[1],vhx[1]});

  plotymin = -1.E-9; plotymax = 3E-9;
  plotyline = x[2];
  Draw({hmcPbdphibkg[2],vhx[2]});

  plotyline = 9999;
  plotymax = 1;

  auto hmcppdphiNSASratio = (TH1F *)hmcppdphiNS->Clone("hmcppdphiNSASratio");
  hmcppdphiNSASratio->Divide(hmcppdphiAS);
  for (int i=0;i<Nb;i++)
    hmcPbdphiNSASratio[i]->Divide(hmcPbdphiNS[i],hmcPbdphiAS[i]);
//    hmcPbdphiNS[i]->Divide(hmcPbdphiAS[i]);


  buildh(Nb+1,0,Nb+1);
  auto hNSASratio = geth("hNSASratio","ratio of NS to AS");
  double e = 0;
  double integral = hmcppdphiNSASratio->IntegralAndError(1, hmcppdphiNSASratio->GetNbinsX(),e);
  cout<<" ? "<<e<<" "<<integral<<endl;
  hNSASratio->SetBinContent(1,hmcppdphiNSASratio->GetBinContent(1));
  hNSASratio->SetBinError(1,hmcppdphiNSASratio->GetBinError(1));
  for (int i=0;i<Nb;i++) {
    e = 0;
    integral = hmcPbdphiNSASratio[i]->IntegralAndError(1, hmcPbdphiNSASratio[i]->GetNbinsX(),e);
    hNSASratio->SetBinContent(i+2,hmcPbdphiNSASratio[i]->GetBinContent(1));
    hNSASratio->SetBinError(i+2,hmcPbdphiNSASratio[i]->GetBinError(1));
  }
  Print(hNSASratio);


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


  cout<<" bkgfractionInNearSide : ";
  for (auto f:bkgfractionInNearSide) {
    cout<<(f < 0 ? 0: f)<<",";
  }
  cout<<endl;

}

void hydjetestimation(bool firstRun = true)
{
  macro m("hydjetestimation_upd",firstRun);
  dphisignalconsistency();
}