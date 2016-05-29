#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

void PrintVector(TString label, vector<float> v)
{
  cout<<" "<<label<<" \t: ";
  if (v.size()==0) {
    cout<<endl;
    return;
  }
  for (unsigned i=0;i<v.size()-1;i++)
    cout<<(v[i] < 0 ? 0: v[i])<<",";
  cout<<v[v.size()-1]<<endl;
}

void dphisignalconsistency()
{

  auto fmcpp = config.getfile_djt("mcppbfa");
  auto fmcPb = config.getfile_djt("mcPbbfa");

  seth(3,0,PI);
  auto hmcppdphi = geth("hmcppdphi","pp;#Delta#phi");
  auto hdtppdphi = geth("hdtppdphi","data pp;#Delta#phi");
  seth(1,0,PI);
  auto hmcppdphiNS = geth("hmcppdphiNS","pp NS;#Delta#phi");
  auto hmcppdphiAS = geth("hmcppdphiAS","pp AS;#Delta#phi");

  auto hdtppdphiNS = geth("hdtppdphiNS","data pp NS;#Delta#phi");
  auto hdtppdphiAS = geth("hdtppdphiAS","data pp AS;#Delta#phi");

  // auto hmcppdphiNSASratio = geth("hmcppdphiNSASratio","pp;#Delta#phi;CR/AS");

  int Nb = bins.size()-1;

  vector<TH1 *> hmcPbdphi(Nb),hmcPbdphiNS(Nb),hmcPbdphiAS(Nb),hmcPbdphiNSall(Nb),hmcPbdphiASall(Nb),
                hmcPbdphibkg(Nb),hmcPbdphiNSbkg(Nb),hmcPbdphiASbkg(Nb),hmcPbdphiNSASratio(Nb),hmcPbdphiNSASallratio(Nb),
                hdtPbdphiNS(Nb),hdtPbdphiAS(Nb),hdtPbdphiNSASratio(Nb);
                
  for (int i=0;i<Nb;i++) {
    seth(3,0,PI);
    hmcPbdphi[i] = geth(Form("hmcPbdphi%d",i),Form("Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphibkg[i] = geth(Form("hmcPbdphibkg%d",i),Form("Hydjet Pb %s;#Delta#phi",binnames[i].Data()));
    seth(1,0,PI);
    hmcPbdphiNS[i] = geth(Form("hmcPbdphiNS%d",i),Form("NS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSbkg[i] = geth(Form("hmcPbdphiNSbkg%d",i),Form("NS Background Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSall[i] = geth(Form("hmcPbdphiNSall%d",i),Form("NS ALL Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiAS[i] = geth(Form("hmcPbdphiAS%d",i),Form("AS Signal Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiASbkg[i] = geth(Form("hmcPbdphiASbkg%d",i),Form("AS Background Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiASall[i] = geth(Form("hmcPbdphiASall%d",i),Form("AS ALL Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSASratio[i] = geth(Form("hmcPbdphiNSASratio%d",i),Form("NS AS ratio Pb %s;#Delta#phi",binnames[i].Data()));
    hmcPbdphiNSASallratio[i] = geth(Form("hmcPbdphiNSASallratio%d",i),Form("PbPb %s;#Delta#phi;CR/AS",binnames[i].Data()));

    hdtPbdphiNS[i] = geth(Form("hdtPbdphiNS%d",i),Form("NS Pb data%s;#Delta#phi",binnames[i].Data()));
    hdtPbdphiAS[i] = geth(Form("hdtPbdphiAS%d",i),Form("AS Pb data%s;#Delta#phi",binnames[i].Data()));
    hdtPbdphiNSASratio[i] = geth(Form("hdtPbdphiNSASratio%d",i),Form("PbPb %s;#Delta#phi;CR/AS",binnames[i].Data()));
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


  // DATA DRIVEN ESTIMATION
  auto fdtpp = config.getfile_djt("dtppjpf");
  Fill(fdtpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL"},[&] (dict d) {
    float w = d["weight"];

    if (d["jtpt1"]>pt1cut && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>pt2cut) {//&& d["discr_csvV1_1"]>0.9
      hdtppdphi->Fill(d["dphiSL1"],w);
      if (NearSide(d)) hdtppdphiNS->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hdtppdphiAS->Fill(d["dphiSL1"],w);
    }

  });

  auto fdtPb = config.getfile_djt("dtPbbjt");
  Fill(fdtPb,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","bin","discr_csvV1_1","jteta1","jtetaSL"},[&] (dict d) {
    float w = d["weight"];

    int b = getbinindex(d["bin"]);
    
    if (d["jtpt1"]>pt1cut && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>pt2cut) { //d["discr_csvV1_1"]>0.9
      if (NearSide(d)) hdtPbdphiNS[b]->Fill(d["dphiSL1"],w);
      if (AwaySide(d)) hdtPbdphiAS[b]->Fill(d["dphiSL1"],w);
    }

  });

  vector<float> truex, centralx, datax;

  vector<float> x;//, y;
  vector<TH1F *>vhx;

  vector<float> NStoASall;

  float PbPb1sig = hmcPbdphiNS[2]->GetBinContent(1);//estimate from central signal
  float PbPb3sig = hmcPbdphiAS[2]->GetBinContent(1);//estimate from central signal

  // float dataNStoAS = hmcppdphiAS->GetBinContent(1)/hmcppdphiNS->GetBinContent(1);

  float pp1 = hmcppdphiNS->GetBinContent(1);
  float pp3 = hmcppdphiAS->GetBinContent(1);
  NStoASall.push_back(pp3/pp1);

  float pp1data = hdtppdphiNS->GetBinContent(1);
  float pp3data = hdtppdphiAS->GetBinContent(1);

  // hmcppdphiNSASratio->Divide(hmcppdphiNS[i],hmcppdphiAS[i]);

  for (int i=0;i<Nb;i++) {
    
    float PbPb1 = hmcPbdphiNSall[i]->GetBinContent(1);
    float PbPb3 = hmcPbdphiASall[i]->GetBinContent(1);
    hmcPbdphiNSASallratio[i]->Divide(hmcPbdphiNSall[i],hmcPbdphiASall[i]);

    float PbPb1bkg = hmcPbdphiNSbkg[i]->GetBinContent(1);
    truex.push_back(PbPb1bkg/PbPb1);

    seth(1,0,PI);
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

    hxx->Multiply(hmcPbdphiNSall[i]);
    vhx.push_back(hxx);


    //just for fun... turn on quenching of the away-side by alpha
    double alpha = 0.3;
    double Bprime = 1/alpha*B;
    float xalpha_ = (1-PbPb3/PbPb1*Bprime)/(1-Bprime);
    cout<<" bkg fraction in NS = "<<xx_<<", quenched by "<<alpha<<" bkg = "<<xalpha_<<endl;
    cout<<setprecision(10)<<xx_*PbPb1<<endl;

    x.push_back(xx_*PbPb1);

    bkgfractionInNearSide[i] = xx_;

    //estimate from PbPb signal
    float Bpbpb = PbPb1sig/PbPb3sig;
    float xxpbpb_ = (1-PbPb3/PbPb1*Bpbpb)/(1-Bpbpb);

    centralx.push_back(xxpbpb_);

    float PbPb1data = hdtPbdphiNS[i]->GetBinContent(1);
    float PbPb3data = hdtPbdphiAS[i]->GetBinContent(1);
    hdtPbdphiNSASratio[i]->Divide(hdtPbdphiNS[i],hdtPbdphiAS[i]);

    //estimate from data
    float Bdata = pp1data/pp3data;
    float xxdata_ = (1-PbPb3data/PbPb1data*Bdata)/(1-Bdata);

    datax.push_back(xxdata_);

  }



  plotlegendpos = TopLeft;
  plotymin = 0; plotymax = 1E-7;
  plotyline = x[0];
  // Draw({hmcPbdphibkg[0],vhx[0]});

  plotymin = 0; plotymax = 1E-8;
  plotyline = x[1];
  // Draw({hmcPbdphibkg[1],vhx[1]});

  plotymin = -1.E-9; plotymax = 3E-9;
  plotyline = x[2];
  // Draw({hmcPbdphibkg[2],vhx[2]});

  plotyline = 9999;
  plotymax = 1;

  auto hdtppdphiNSASratio = (TH1F *)hdtppdphiNS->Clone("hdtppdphiNSASratio");
  hdtppdphiNSASratio->Divide(hdtppdphiAS);

  auto hmcppdphiNSASratio = (TH1F *)hmcppdphiNS->Clone("hmcppdphiNSASratio");
  hmcppdphiNSASratio->Divide(hmcppdphiAS);
  for (int i=0;i<Nb;i++)
    hmcPbdphiNSASratio[i]->Divide(hmcPbdphiNS[i],hmcPbdphiAS[i]);
//    hmcPbdphiNS[i]->Divide(hmcPbdphiAS[i]);


  seth(Nb+1,0,Nb+1);
  auto hNSASratio = geth("hNSASratio","ratio of NS to AS;;control region/away-side signal yields ratio");

  auto hbkgfraction = geth("hbkgfraction","MC, estimated from Pythia;;Bkg fraction in CR");
  auto hbkgfractiontrue = geth("hbkgfractiontrue","MC, true fraction;;Bkg fraction in CR");
  auto hbkgfractioncent = geth("hbkgfractioncent","MC, estimated from Pythia+Hydjet;;Bkg fraction in CR");
  auto hbkgfractiondata = geth("hbkgfractiondata","Data, estimated from pp;;Bkg fraction in CR");
  auto hNSASallratio = geth("hNSASallratio","MC;;control region/away-side yields ratio");
  auto hNSASratiodata = geth("hNSASratiodata","Data;;control region/away-side yields ratio");


  double e = 0;
  double integral = hmcppdphiNSASratio->IntegralAndError(1, hmcppdphiNSASratio->GetNbinsX(),e);
  cout<<" ? "<<e<<" "<<integral<<endl;
  hNSASratio->SetBinContent(1,hmcppdphiNSASratio->GetBinContent(1));
  hNSASratio->SetBinError(1,hmcppdphiNSASratio->GetBinError(1));
  hNSASallratio->SetBinContent(1,hmcppdphiNSASratio->GetBinContent(1));
  hNSASallratio->SetBinError(1,hmcppdphiNSASratio->GetBinError(1));

  hNSASratiodata->SetBinContent(1,hdtppdphiNSASratio->GetBinContent(1));
  hNSASratiodata->SetBinError(1,hdtppdphiNSASratio->GetBinError(1));

  hbkgfraction->SetBinContent(1,0);
  hbkgfraction->SetBinError(1,1E-3);
  hbkgfractiontrue->SetBinContent(1,0);
  hbkgfractiontrue->SetBinError(1,1E-3);
  hbkgfractioncent->SetBinContent(1,0);
  hbkgfractioncent->SetBinError(1,1E-3);
  hbkgfractiondata->SetBinContent(1,0);
  hbkgfractiondata->SetBinError(1,1E-3);

  for (int i=0;i<Nb;i++) {
    e = 0;
    integral = hmcPbdphiNSASratio[i]->IntegralAndError(1, hmcPbdphiNSASratio[i]->GetNbinsX(),e);
    hNSASratio->SetBinContent(i+2,hmcPbdphiNSASratio[i]->GetBinContent(1));
    hNSASratio->SetBinError(i+2,hmcPbdphiNSASratio[i]->GetBinError(1));

    hbkgfraction->SetBinContent(i+2,bkgfractionInNearSide[i]);
    hbkgfraction->SetBinError(i+2,1E-3);

    hbkgfractiontrue->SetBinContent(i+2,truex[i]);
    hbkgfractiontrue->SetBinError(i+2,1E-3);    

    hbkgfractioncent->SetBinContent(i+2,centralx[i]);
    hbkgfractioncent->SetBinError(i+2,1E-3);

    hbkgfractiondata->SetBinContent(i+2,datax[i]);
    hbkgfractiondata->SetBinError(i+2,1E-3);
        
    hNSASallratio->SetBinContent(i+2,hmcPbdphiNSASallratio[i]->GetBinContent(1));
    hNSASallratio->SetBinError(i+2,hmcPbdphiNSASallratio[i]->GetBinError(1));
    hNSASratiodata->SetBinContent(i+2,hdtPbdphiNSASratio[i]->GetBinContent(1));
    hNSASratiodata->SetBinError(i+2,hdtPbdphiNSASratio[i]->GetBinError(1));

    
  }


  vector<TString> axnames;
  axnames.push_back("pp");
  for (auto s:binnames) axnames.push_back(s);
  RenameBinLabelsX(hNSASratio,axnames);
  RenameBinLabelsX(hNSASratiodata,axnames);
  RenameBinLabelsX(hbkgfraction,axnames);
  RenameBinLabelsX(hbkgfractiontrue,axnames);  
  RenameBinLabelsX(hbkgfractioncent,axnames);
  RenameBinLabelsX(hbkgfractiondata,axnames);
  RenameBinLabelsX(hNSASallratio,axnames);

  hmcPbdphi.push_back(hmcppdphi);
  Normalize(hmcPbdphi);

  plotlegendpos = TopLeft;
  // Draw(hmcPbdphi);


  ShuffleBins(hNSASratio,{1,4,3,2});
  ShuffleBins(hbkgfraction,{1,4,3,2});
  ShuffleBins(hbkgfractiontrue,{1,4,3,2});  
  ShuffleBins(hbkgfractioncent,{1,4,3,2});  
  ShuffleBins(hbkgfractiondata,{1,4,3,2});  
  ShuffleBins(hNSASallratio,{1,4,3,2}); 
  ShuffleBins(hNSASratiodata,{1,4,3,2});


  plotymin = 0;
  plotymax = 0.15;
  plotlegendpos = None;
  plotyline = hNSASratio->GetBinContent(1);
  Draw({hNSASratio});
  plotyline = 9999;

  plotymax = 0.5;
    plotlegendpos = TopLeft;
  SetMC({hNSASallratio});
  Draw({hNSASratiodata,hNSASallratio});

  plotymax = 1;
  plotlegendpos = TopLeft;
  hbkgfractiontrue->SetMarkerStyle(kOpenSquare);
  Draw({hbkgfractiontrue,hbkgfraction,hbkgfractioncent},"E1");
  Draw({hbkgfractiontrue,hbkgfraction,hbkgfractiondata},"E1");

  PrintVector("true fraction in NS",truex);

  PrintVector("bkgfractionInNearSide",bkgfractionInNearSide);


}

void hydjetestimation(bool firstRun = true)
{
  macro m("hydjetestimation",firstRun);
  dphisignalconsistency();
}