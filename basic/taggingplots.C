#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

int getind(int refparton, bool isSignal)
{
  if (isSignal) {
    if (abs(refparton)==5) return 0; else
    if (abs(refparton)==4) return 1; else
    return 3;
  } 
  else return 2;
}

vector<TGraph *> roc;
vector<TGraph *> roccsvcut;


void makeplots(int binMin, int binMax, float fr)
{
  bool pp = binMin==-1 && binMax==-1; //pp = false, PbPb = true

  auto fmcqcd = config.getfile_djt(pp ? "mcppqcd":"mcPbqcd");
  auto fmcbfa = config.getfile_djt(pp ? "mcppbfa":"mcPbbfa");
  auto fdtinc = config.getfile_djt(pp ? "dtppjpf":"dtPbjcl");
  auto fdtbjt = config.getfile_djt("dtPbbjt");

  setv({"b","c","UE","usdg"});
  vector<int> colors = {darkred, darkgreen, kYellow,darkblue};

  buildNamesuffix = pp? "pp" : Form("_bin%d%d",binMin, binMax);
  seth(20,-0.5,19.5); buildndiv = 20;
  auto hmcnumtag = geth("hmcnumtag","MC;Number of tagged jets");
  auto hdtnumtag = geth("hdtnumtag","Data;Number of tagged jets");
  buildndiv = -1;

  seth(10,0,1);
  auto hmcLJdiscr = getv("hmcLJdiscr","MC;Leading jet CSV");
  auto hdtLJdiscr = geth("hdtLJdiscr","Data;Leading jet CSV");
  auto hmcPJdiscr = getv("hmcPJdiscr","MC;Subleading jet CSV");
  auto hdtPJdiscr = geth("hdtPJdiscr","Data;Subleading jet CSV");

  auto hmcLJdiscr_ntrkGT2 = getv("hmcLJdiscr_ntrkGT2","MC;Leading jet CSV");
  auto hdtLJdiscr_ntrkGT2 = geth("hdtLJdiscr_ntrkGT2","Data;Leading jet CSV");
  auto hmcPJdiscr_ntrkGT2 = getv("hmcPJdiscr_ntrkGT2","MC;Subleading jet CSV");
  auto hdtPJdiscr_ntrkGT2 = geth("hdtPJdiscr_ntrkGT2","Data;Subleading jet CSV");


  seth(5,-0.5,4.5); buildndiv = 5;
  auto hmcLJnsvtx = getv("hmcLJnsvtx","MC;Leading jet # of SV");
  auto hdtLJnsvtx = geth("hdtLJnsvtx","Data;Leading jet # of SV");
  auto hmcPJnsvtx = getv("hmcPJnsvtx","MC;Subleading jet # of SV");
  auto hdtPJnsvtx = geth("hdtPJnsvtx","Data;Subleading jet # of SV");
  buildndiv = -1;


  seth(20,0,80);
  auto hmcLJsvtxdls = getv("hmcLJsvtxdls","MC;Leading jet SV distance significance");
  auto hdtLJsvtxdls = geth("hdtLJsvtxdls","Data;Leading jet SV distance significance");
  auto hmcPJsvtxdls = getv("hmcPJsvtxdls","MC;Subleading jet SV distance significance");
  auto hdtPJsvtxdls = geth("hdtPJsvtxdls","Data;Subleading jet SV distance significance");

  seth(12,0,6);
  auto hmcLJsvtxm = getv("hmcLJsvtxm","MC;Leading jet SV mass [GeV]");
  auto hdtLJsvtxm = geth("hdtLJsvtxm","Data;Leading jet SV mass [GeV]");
  auto hmcPJsvtxm = getv("hmcPJsvtxm","MC;Subleading jet SV mass [GeV]");
  auto hdtPJsvtxm = geth("hdtPJsvtxm","Data;Subleading jet SV mass [GeV]");

  auto hmcLJsvtxm_ntrkGT2 = getv("hmcLJsvtxm_ntrkGT2","MC;Leading jet SV mass [GeV]");
  auto hdtLJsvtxm_ntrkGT2 = geth("hdtLJsvtxm_ntrkGT2","Data;Leading jet SV mass [GeV]");
  auto hmcPJsvtxm_ntrkGT2 = getv("hmcPJsvtxm_ntrkGT2","MC;Subleading jet SV mass [GeV]");
  auto hdtPJsvtxm_ntrkGT2 = geth("hdtPJsvtxm_ntrkGT2","Data;Subleading jet SV mass [GeV]");

  seth(10,0,100);//40
  auto hmcLJsvtxpt = getv("hmcLJsvtxpt","MC;Leading jet SV p_{T} [GeV]");
  auto hdtLJsvtxpt = geth("hdtLJsvtxpt","Data;Leading jet SV p_{T} [GeV]");
  auto hmcPJsvtxpt = getv("hmcPJsvtxpt","MC;Subleading jet SV p_{T} [GeV]");
  auto hdtPJsvtxpt = geth("hdtPJsvtxpt","Data;Subleading jet SV p_{T} [GeV]");

  seth(100,-0.1,0.1);//40
  auto hmcLJip3d = getv("hmcLJip3d","MC;Leading jet SV p_{T} [GeV]");
  auto hdtLJip3d = geth("hdtLJip3d","Data;Leading jet SV p_{T} [GeV]");
  auto hmcPJip3d = getv("hmcPJip3d","MC;Subleading jet SV p_{T} [GeV]");
  auto hdtPJip3d = geth("hdtPJip3d","Data;Subleading jet SV p_{T} [GeV]");


  vector<float> csvrange = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95};
  int ncsv = csvrange.size();

  seth(5,0,5);
  vector<vector<TH1F *> > vhPairCode (ncsv); //for each csv1, csv2
  vector<float> effnum(ncsv);
  float effden = 0;

  for (int i=0;i<ncsv;i++)
    for (int j = 0;j<ncsv;j++)
      vhPairCode[i].push_back(geth(Form("vhPairCode_%.2f_%.2f",csvrange[i],csvrange[j])));

  auto fillmchistsfunc = [&] (dict &d) {
    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;
    if (d["pthat"]<pthatcut) return;

    float w = weight12(d);

    hmcnumtag->Fill(d["numTagged"],w);
    if (d["numTagged"]>6) return;

    int LJind = getind(d["refparton_flavorForB1"],true);
    // int PJind = getind(d["refparton_flavorForBSL"],IsSignal(d));
    int PJind = getind(d["refparton_flavorForB2"],d["subid2"]==0 && d["refpt2"]>20);//IsSignal(d));

    //Leading jet

    if (!(d["jtpt1"]>pt1cut && d["refpt1"]>50)) return;

    hmcLJdiscr[LJind]->Fill(d["discr_csvV1_1"],w);
    hmcLJnsvtx[LJind]->Fill(d["nsvtx1"],w);
    hmcLJsvtxdls[LJind]->Fill(d["svtxdls1"],w);
    hmcLJsvtxm[LJind]->Fill(d["svtxm1"],w);
    hmcLJsvtxpt[LJind]->Fill(d["svtxpt1"],w);
    if (d["svtxntrk1"]>1){
      hmcLJsvtxm_ntrkGT2[LJind]->Fill(d["svtxm1"],w);
      hmcLJdiscr_ntrkGT2[LJind]->Fill(d["discr_csvV1_1"],w);
    }

    //partner jet

    if (d["discr_csvV1_1"]>0.9 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hmcPJdiscr[PJind]->Fill(d["discr_csvV1_2"],w);
      hmcPJnsvtx[PJind]->Fill(d["nsvtx2"],w);
      hmcPJsvtxdls[PJind]->Fill(d["svtxdls2"],w);
      hmcPJsvtxm[PJind]->Fill(d["svtxm2"],w);
      hmcPJsvtxpt[PJind]->Fill(d["svtxpt2"],w);
      if (d["svtxntrk2"]>1) {
        hmcPJsvtxm_ntrkGT2[PJind]->Fill(d["svtxm2"],w);
        hmcPJdiscr_ntrkGT2[PJind]->Fill(d["discr_csvV1_2"],w);
      }
    }

    if (d["jtptSignal2"]>pt2cut && d["dphiSignal21"]>PI23) {
      float ws2 = weight1Signal2(d);

      float csv1 = d["discr_csvV1_1"];
      float csv2 = d["discr_csvV1_Signal2"];
      float paircode = d["pairCodeSignal21"];//????????
  
      for (int i=0;i<ncsv;i++)
        for (int j = 0;j<ncsv;j++)
          if (csv1>csvrange[i] && csv2>csvrange[j])
            vhPairCode[i][j]->Fill(paircode,ws2);
  

      if (paircode==0) {
        effden+=ws2;
        for (int i=0;i<ncsv;i++)
          if (csv1>csvrange[i] && csv2>csvrange[i])
            effnum[i]+=ws2;
      }


    }

  };


  Fill(fmcbfa,[&] (dict &d) {

    if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForBSL"])!=5) return; //use bjt sample with b on either side 
    // if (abs(d["refparton_flavorForB1"])!=5) return; //use bjt sample with b on either side 
    if (d["event"]==128751 || d["event"]==1551232 || d["event"]==3043866) return;

    fillmchistsfunc(d);

  },fr);

  Fill(fmcqcd,[&] (dict &d) {

    if (abs(d["refparton_flavorForB1"])==5 || abs(d["refparton_flavorForBSL"])==5) return; //use qcd sample with no b on either side 
    // if (abs(d["refparton_flavorForB1"])==5) return; //use qcd sample with no b on either side 
    if (d["event"]==128751) return;


    fillmchistsfunc(d);

  },fr);

  Fill(fdtinc,[&] (dict &d) {

    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];

    hdtnumtag->Fill(d["numTagged"],w);
    if (d["numTagged"]>6) return;

    if (d["jtpt1"]>pt1cut){
      hdtLJdiscr->Fill(d["discr_csvV1_1"],w);
      hdtLJnsvtx->Fill(d["nsvtx1"],w);
      hdtLJsvtxdls->Fill(d["svtxdls1"],w);
      hdtLJsvtxm->Fill(d["svtxm1"],w);
      hdtLJsvtxpt->Fill(d["svtxpt1"],w);
      if (d["svtxntrk1"]>1){
        hdtLJsvtxm_ntrkGT2->Fill(d["svtxm1"],w);
        hdtLJdiscr_ntrkGT2->Fill(d["discr_csvV1_1"],w);
      }
    }

    // if (d["jtpt1"]>pt1cut && d[jtptSL]>pt2cut && d[dphiSL1]>PI23){
    //   hdtPJdiscr->Fill(d["discr_csvV1_SL"],w);
    //   hdtPJnsvtx->Fill(d["nsvtxSL"],w);
    //   hdtPJsvtxdls->Fill(d["svtxdlsSL"],w);
    //   hdtPJsvtxm->Fill(d["svtxmSL"],w);
    //   if (d["svtxntrkSL"]>1) {
    //     hdtPJsvtxm_ntrkGT2->Fill(d["svtxmSL"],w);
    //     hdtPJdiscr_ntrkGT2->Fill(d["discr_csvV1_SL"],w);
    //   }
    // }

    // fill subleading jet in PbPb with b-jet sample only
    if (pp && d["jtpt1"]>pt1cut && d["discr_csvV1_1"]>0.9 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hdtPJdiscr->Fill(d["discr_csvV1_2"],w);
      hdtPJnsvtx->Fill(d["nsvtx2"],w);
      hdtPJsvtxdls->Fill(d["svtxdls2"],w);
      hdtPJsvtxm->Fill(d["svtxm2"],w);
      hdtPJsvtxpt->Fill(d["svtxpt2"],w);
      if (d["svtxntrk2"]>1) {
        hdtPJsvtxm_ntrkGT2->Fill(d["svtxm2"],w);
        hdtPJdiscr_ntrkGT2->Fill(d["discr_csvV1_2"],w);
      }
    }

  },fr);

  if (!pp) //in PbPb case only use specific b-jet sample
  Fill(fdtbjt,[&] (dict &d) {

    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];

    // hdtnumtag->Fill(d["numTagged"],w);
    // if (d["numTagged"]>6) return;

    if (d["jtpt1"]>pt1cut && d["discr_csvV1_1"]>0.9 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hdtPJdiscr->Fill(d["discr_csvV1_2"],w);
      hdtPJnsvtx->Fill(d["nsvtx2"],w);
      hdtPJsvtxdls->Fill(d["svtxdls2"],w);
      hdtPJsvtxm->Fill(d["svtxm2"],w);
      hdtPJsvtxpt->Fill(d["svtxpt2"],w);
      if (d["svtxntrk2"]>1) {
        hdtPJsvtxm_ntrkGT2->Fill(d["svtxm2"],w);
        hdtPJdiscr_ntrkGT2->Fill(d["discr_csvV1_2"],w);
      }
    }

  },fr);

  MakeOverflowVisibleAll();

  auto hsmcLJdiscr = stackhists(hmcLJdiscr, colors, "hsmcLJdiscr");
  auto hsmcPJdiscr = stackhists(hmcPJdiscr, colors, "hsmcPJdiscr");
  auto hsmcLJdiscr_ntrkGT2 = stackhists(hmcLJdiscr_ntrkGT2, colors, "hsmcLJdiscr_ntrkGT2");
  auto hsmcPJdiscr_ntrkGT2 = stackhists(hmcPJdiscr_ntrkGT2, colors, "hsmcPJdiscr_ntrkGT2");
  auto hsmcLJnsvtx = stackhists(hmcLJnsvtx, colors, "hsmcLJnsvtx");
  auto hsmcPJnsvtx = stackhists(hmcPJnsvtx, colors, "hsmcPJnsvtx");
  auto hsmcLJnsvtpt = stackhists(hmcLJsvtxpt, colors, "hsmcLJnsvtpt");
  auto hsmcPJnsvtpt = stackhists(hmcPJsvtxpt, colors, "hsmcPJnsvtpt");
  auto hsmcLJsvtxdls = stackhists(hmcLJsvtxdls, colors, "hsmcLJsvtxdls");
  auto hsmcPJsvtxdls = stackhists(hmcPJsvtxdls, colors, "hsmcPJsvtxdls");
  auto hsmcLJsvtxm = stackhists(hmcLJsvtxm, colors, "hsmcLJsvtxm");
  auto hsmcPJsvtxm = stackhists(hmcPJsvtxm, colors, "hsmcPJsvtxm");
  auto hsmcLJsvtxpt = stackhists(hmcLJsvtxpt, colors, "hsmcLJsvtxpt");
  auto hsmcPJsvtxpt = stackhists(hmcPJsvtxpt, colors, "hsmcPJsvtxpt");
  auto hsmcLJsvtxm_ntrkGT2 = stackhists(hmcLJsvtxm_ntrkGT2, colors, "hsmcLJsvtxm_ntrkGT2");
  auto hsmcPJsvtxm_ntrkGT2 = stackhists(hmcPJsvtxm_ntrkGT2, colors, "hsmcPJsvtxm_ntrkGT2");



  SetData({hdtLJdiscr_ntrkGT2,hdtPJdiscr_ntrkGT2,hdtLJdiscr,hdtPJdiscr,hdtLJnsvtx,hdtPJnsvtx,hdtnumtag,hdtLJsvtxdls,hdtPJsvtxdls,hdtLJsvtxm,hdtPJsvtxm,hdtLJsvtxpt,hdtPJsvtxpt,hdtLJsvtxm_ntrkGT2,hdtPJsvtxm_ntrkGT2});
  SetB({hdtLJdiscr_ntrkGT2,hdtPJdiscr_ntrkGT2,hdtLJdiscr,hdtPJdiscr,hdtLJnsvtx,hdtPJnsvtx,hdtnumtag,hdtLJsvtxdls,hdtPJsvtxdls,hdtLJsvtxm,hdtPJsvtxm,hdtLJsvtxpt,hdtPJsvtxpt,hdtLJsvtxm_ntrkGT2,hdtPJsvtxm_ntrkGT2});
  
  SetMC({hmcnumtag});
  SetB({hmcnumtag});

  plotylog = true;
  plotytitle = "Event fractions";
  // aktstring = "anti-k_{T} Pu R=0.4 |#eta|<2.0";  
  aktstring = pp ? "pp" : Form("PbPb %d-%d %%",binMin/2,binMax/2);

  plotlegendorder = {3,2,0};
  plotlegenddx = 0.1;
  textposx = 0.21;

  plotytitle = "Counts";

  Normalize({hdtnumtag,hmcnumtag});
  DrawCompare(hdtnumtag,hmcnumtag);

  plotsecondline = Form("p_{T,1}>%d GeV", (int)pt1cut);
  plotthirdline = "";
  //plotymin = 1;
  DrawCompare(hdtLJdiscr,hsmcLJdiscr);
  plotymin = 9999;
  plotymax  = pp? 1E8:1E7;
  DrawCompare(hdtLJnsvtx,hsmcLJnsvtx);
  plotymax  = 9999;
  DrawCompare(hdtLJsvtxdls,hsmcLJsvtxdls);
  // plotylog = false;
  DrawCompare(hdtLJsvtxm,hsmcLJsvtxm);
  plotylog = true;
  DrawCompare(hdtLJsvtxpt,hsmcLJsvtxpt);
  plotylog = false;

  // plotylog = true;
  // plotthirdline = "# s.v. tracks>2";
  // DrawCompare(hdtLJdiscr_ntrkGT2,hsmcLJdiscr_ntrkGT2);
  // plotylog = false;
  // DrawCompare(hdtLJsvtxm_ntrkGT2,hsmcLJsvtxm_ntrkGT2);
  // plotylog = true;

  // aktstring = "anti-k_{T} Pu R=0.4 |#eta|<2.0";
  if (!pp)
    plotlegendorder = {3,2,0,1};
  
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";
  // plotymin = 0.001;
  //plotymin = 1;
    plotylog = true;
  DrawCompare(hdtPJdiscr,hsmcPJdiscr);
  plotymin = 9999;
  plotymax  = pp? 1E8:1E6;
  DrawCompare(hdtPJnsvtx,hsmcPJnsvtx);
  plotymax  =9999;
  DrawCompare(hdtPJsvtxdls,hsmcPJsvtxdls);
  // plotylog = false;
  DrawCompare(hdtPJsvtxm,hsmcPJsvtxm);

  plotylog = true;
  DrawCompare(hdtPJsvtxpt,hsmcPJsvtxpt);
  plotylog = false;

  // plotylog = true;
  // plotthirdline += "; # s.v. tracks>2";
  // DrawCompare(hdtPJdiscr_ntrkGT2,hsmcPJdiscr_ntrkGT2);
  // plotylog = false;
  // DrawCompare(hdtPJsvtxm_ntrkGT2,hsmcPJsvtxm_ntrkGT2);


  vector<vector<float> > pur(ncsv);

  for (int i=0;i<ncsv;i++)
    for (int j=0;j<ncsv;j++)
      pur[i].push_back(vhPairCode[i][j]->GetBinContent(1)/vhPairCode[i][j]->Integral());

  float effcsvcut = 0,purcsvcut = 0;
  vector<float> eff12,pur12;
  for (int i=0;i<ncsv;i++) {
    pur12.push_back(pur[i][i]);
    eff12.push_back(effnum[i]/effden);

    if (csvrange[i]==csvcut1) {
      effcsvcut = eff12[i];
      purcsvcut = pur12[i];
    }

  }

  // pur12.push_back(1);
  // eff12.push_back(0);


  vector<TGraph *> gg;
  auto c = getc();
  c->SetGrid();
  plotlegendpos = TopLeft;
  auto l = new TLegend(0.18,0.4,0.4,0.8);
  for (int i=0;i<ncsv;i++) {
    auto g = new TGraph(ncsv,&csvrange[0],&(pur[i][0]));
    g->SetMinimum(0);
    g->SetMaximum(1);
    g->SetMarkerColor(TColor::GetColorDark(i>=10 ? i+1 : i));
    gg.push_back(g);
    l->AddEntry(g,Form("csv1>%.2f",csvrange[i]),"P");
    g->Draw(i==0 ? "APL" : "PL,same");
    g->GetXaxis()->SetTitle("csv2");
  }
  l->Draw();
  SavePlot(c,binMin==-1 ? "purity_pp" : Form("purity_%d_%d",binMin,binMax));



  // auto croc = getc();
  // croc->SetGrid();
  auto groc = new TGraph(ncsv,&pur12[0],&eff12[0]);
  groc->SetTitle(pp ? "pp" : Form("%d-%d%%",(int)binMin/2,(int)binMax/2));
  // groc->Draw("APL");
  groc->GetXaxis()->SetTitle("b-tagging purity");
  groc->GetYaxis()->SetTitle("b-tagging efficiency");
  // plotlegendpos = TopLeft;
  // SavePlot(croc,binMin==-1 ? "roc_pp" : Form("roc_%d_%d",binMin,binMax));

  roc.push_back(groc);

  auto gcsvroc = new TGraph(1,&purcsvcut, &effcsvcut);
  gcsvroc->SetMarkerColor(kRed);
  gcsvroc->SetMarkerSize(1.3);
  roccsvcut.push_back(gcsvroc);

}

void WriteROC()
{

  auto f = new TFile(plotfoldername+"/ROC.root","recreate");
  for (unsigned i=0;i<roc.size();i++) {
    roc[i]->Write(Form("roc%d",i));
    roccsvcut[i]->Write(Form("roccsvcut%d",i));
  }
  f->Close();

}

void ReadROC()
{

  auto f = new TFile(plotfoldername+"/ROC.root");
  for (unsigned i=0;i<4;i++) {
    roc.push_back((TGraph *)f->Get(Form("roc%d",i)));
    roccsvcut.push_back((TGraph *)f->Get(Form("roccsvcut%d",i)));
  }
  f->Close();
  cout<<"ROC curves read!"<<endl;

}

void taggingplots()
{
  macro m("taggingplots_0824",true);
  
  float fr = 1.;

  makeplots(-1,-1,fr);
  makeplots(0,20,fr);

  // makeplots(20,60,fr);
  // makeplots(60,200,fr);


//read ROC from file
 ReadROC();

  auto c = getc();
  // c->SetGrid();
  plotlegenddx = 0.03;;
  plotlegendpos = TopRight;
  auto l = getLegend();

roc[0]->SetTitle("pp");
roc[3]->SetTitle("PbPb 30-100%");
roc[2]->SetTitle("PbPb 10-30%");
roc[1]->SetTitle("PbPb 0-10%");

  for (unsigned i=0;i<roc.size();i++) {
    roc[i]->SetMarkerColor(TColor::GetColorDark(i+2));
    roc[i]->SetLineColor(TColor::GetColorDark(i+2));
    roc[i]->SetLineStyle(7);
    roc[i]->SetMarkerStyle(i==3 ? 32 : 24+i);    
    roc[i]->Draw(i==0 ? "APL" : "PL,same");
    // l->AddEntry(roc[i],roc[i]->GetTitle(),"P");

    roc[i]->SetMinimum(0.);
    roc[i]->SetMaximum(1);

    roc[i]->GetXaxis()->SetTitle("b-tagging purity");
    roc[i]->GetYaxis()->SetTitle("b-tagging efficiency");
roc[i]->GetXaxis()->CenterTitle();
roc[i]->GetYaxis()->CenterTitle();

    roc[i]->GetXaxis()->SetRangeUser(0,1);

  roccsvcut[i]->SetMarkerColor(TColor::GetColorDark(i+2));
  roccsvcut[i]->SetMarkerSize(1.5);//2.3);
  roccsvcut[i]->SetMarkerStyle(20+i);  //(29);
  }
  for (unsigned i=0;i<roc.size();i++)
    roccsvcut[i]->Draw("P,same");

  // for (auto g:roc) {
  //   // auto sp = new TSpline3(TString("grs")+g->GetName(),g,"",10,10);
  //   // sp->SetLineColor(kRed);
  //   // sp->Draw("same");

  // auto gs = new TGraphSmooth("normal");
  // auto grout = gs->Approx(g,"linear");
  // grout->Draw("same");

  // }


  l->AddEntry(roc[0],roc[0]->GetTitle(),"P");
  l->AddEntry(roc[3],roc[3]->GetTitle(),"P");
  l->AddEntry(roc[2],roc[2]->GetTitle(),"P");
  l->AddEntry(roc[1],roc[1]->GetTitle(),"P");
  l->SetHeader("Pythia6 (+Hydjet)");

  l->Draw();



  auto ln = new TLine(0,1,1,0);
  ln->SetLineColor(kGray);
  ln->SetLineWidth(1);
  ln->SetLineStyle(7);
  // ln->Draw();

  lumi_sqrtS="";
  extraText = "Simulation";
  CMS_lumi(c, iPeriod, iPos ); 

  SavePlot(c,"ROC");

  // WriteROC();

}



