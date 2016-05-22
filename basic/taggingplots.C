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

void makeplots(int binMin, int binMax, float fr)
{
  bool pp = binMin==-1 && binMax==-1; //pp = false, PbPb = true

  auto fmcqcd = config.getfile_djt(pp ? "mcppqcd":"mcPbqcd");
  auto fmcbfa = config.getfile_djt(pp ? "mcppbfa":"mcPbbfa");
  auto fdtinc = config.getfile_djt(pp ? "dtppjpf":"dtPbj60");
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


  auto fillmchistsfunc = [&] (dict &d) {
    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;
    if (d["pthat"]<pthatcut) return;

    float w = pp ? d["weight"] : weight1SLPbPb(d);

    hmcnumtag->Fill(d["numTagged"],w);
    if (d["numTagged"]>6) return;

    int LJind = getind(d["refparton_flavorForB1"],true);
    // int PJind = getind(d["refparton_flavorForBSL"],IsSignal(d));
    int PJind = getind(d["refparton_flavorForB2"],d["subid2"]==0 && d["refpt2"]>20);//IsSignal(d));


    if (d["jtpt1"]>pt1cut && d["refpt1"]>50){
      hmcLJdiscr[LJind]->Fill(d["discr_csvV1_1"],w);
      hmcLJnsvtx[LJind]->Fill(d["nsvtx1"],w);
      hmcLJsvtxdls[LJind]->Fill(d["svtxdls1"],w);
      hmcLJsvtxm[LJind]->Fill(d["svtxm1"],w);
      if (d["svtxntrk1"]>1){
        hmcLJsvtxm_ntrkGT2[LJind]->Fill(d["svtxm1"],w);
        hmcLJdiscr_ntrkGT2[LJind]->Fill(d["discr_csvV1_1"],w);
      }
    }

    if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && d["discr_csvV1_1"]>0.9 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hmcPJdiscr[PJind]->Fill(d["discr_csvV1_2"],w);
      hmcPJnsvtx[PJind]->Fill(d["nsvtx2"],w);
      hmcPJsvtxdls[PJind]->Fill(d["svtxdls2"],w);
      hmcPJsvtxm[PJind]->Fill(d["svtxm2"],w);
      if (d["svtxntrk2"]>1) {
        hmcPJsvtxm_ntrkGT2[PJind]->Fill(d["svtxm2"],w);
        hmcPJdiscr_ntrkGT2[PJind]->Fill(d["discr_csvV1_2"],w);
      }
    }
  };


  Fill(fmcbfa,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL","subidSL","refptSL","pairCodeSL1",
                  "refparton_flavorForBSL","refparton_flavorForB1","numTagged","refparton_flavorForB2","subid2","refpt2",
                  "jtpt2","dphi21","discr_csvV1_2","nsvtx2","svtxdls2","svtxm2","svtxntrk2","svtxm2","discr_csvV1_2",
                   "nsvtx1","svtxdls1","svtxm1","nsvtxSL","svtxdlsSL","svtxmSL","svtxntrk1","discr_csvV1_SL"},[&] (dict &d) {

    // if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForBSL"])!=5) return; //use bjt sample with b on either side 
    if (abs(d["refparton_flavorForB1"])!=5) return; //use bjt sample with b on either side 

    fillmchistsfunc(d);

  },fr);

  Fill(fmcqcd,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL","subidSL","refptSL","pairCodeSL1",
                  "refparton_flavorForBSL","refparton_flavorForB1","numTagged","refparton_flavorForB2","subid2","refpt2",
                  "jtpt2","dphi21","discr_csvV1_2","nsvtx2","svtxdls2","svtxm2","svtxntrk2","svtxm2","discr_csvV1_2",
                   "nsvtx1","svtxdls1","svtxm1","nsvtxSL","svtxdlsSL","svtxmSL","svtxntrk1","discr_csvV1_SL"},[&] (dict &d) {

    // if (abs(d["refparton_flavorForB1"])==5 || abs(d["refparton_flavorForBSL"])==5) return; //use qcd sample with no b on either side 
    if (abs(d["refparton_flavorForB1"])==5) return; //use qcd sample with no b on either side 

    fillmchistsfunc(d);

  },fr);

  Fill(fdtinc,{"weight","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL",
                  "numTagged",
                  "jtpt2","dphi21","discr_csvV1_2","nsvtx2","svtxdls2","svtxm2","svtxntrk2","svtxm2","discr_csvV1_2",
                   "nsvtx1","svtxdls1","svtxm1","nsvtxSL","svtxdlsSL","svtxmSL","svtxntrk1","discr_csvV1_SL"},[&] (dict &d) {

    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];

    hdtnumtag->Fill(d["numTagged"],w);
    if (d["numTagged"]>6) return;

    if (d["jtpt1"]>pt1cut){
      hdtLJdiscr->Fill(d["discr_csvV1_1"],w);
      hdtLJnsvtx->Fill(d["nsvtx1"],w);
      hdtLJsvtxdls->Fill(d["svtxdls1"],w);
      hdtLJsvtxm->Fill(d["svtxm1"],w);
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
      if (d["svtxntrk2"]>1) {
        hdtPJsvtxm_ntrkGT2->Fill(d["svtxm2"],w);
        hdtPJdiscr_ntrkGT2->Fill(d["discr_csvV1_2"],w);
      }
    }

  },fr);

  if (!pp) //in PbPb case only use specific b-jet sample
  Fill(fdtbjt,{"weight","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL",
                  "numTagged",
                  "jtpt2","dphi21","discr_csvV1_2","nsvtx2","svtxdls2","svtxm2","svtxntrk2","svtxm2","discr_csvV1_2",
                   "nsvtx1","svtxdls1","svtxm1","nsvtxSL","svtxdlsSL","svtxmSL","svtxntrk1","discr_csvV1_SL"},[&] (dict &d) {

    if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];

    // hdtnumtag->Fill(d["numTagged"],w);
    // if (d["numTagged"]>6) return;

    if (d["jtpt1"]>pt1cut && d["discr_csvV1_1"]>0.9 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
      hdtPJdiscr->Fill(d["discr_csvV1_2"],w);
      hdtPJnsvtx->Fill(d["nsvtx2"],w);
      hdtPJsvtxdls->Fill(d["svtxdls2"],w);
      hdtPJsvtxm->Fill(d["svtxm2"],w);
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
  auto hsmcLJsvtxdls = stackhists(hmcLJsvtxdls, colors, "hsmcLJsvtxdls");
  auto hsmcPJsvtxdls = stackhists(hmcPJsvtxdls, colors, "hsmcPJsvtxdls");
  auto hsmcLJsvtxm = stackhists(hmcLJsvtxm, colors, "hsmcLJsvtxm");
  auto hsmcPJsvtxm = stackhists(hmcPJsvtxm, colors, "hsmcPJsvtxm");
  auto hsmcLJsvtxm_ntrkGT2 = stackhists(hmcLJsvtxm_ntrkGT2, colors, "hsmcLJsvtxm_ntrkGT2");
  auto hsmcPJsvtxm_ntrkGT2 = stackhists(hmcPJsvtxm_ntrkGT2, colors, "hsmcPJsvtxm_ntrkGT2");



  SetData({hdtLJdiscr_ntrkGT2,hdtPJdiscr_ntrkGT2,hdtLJdiscr,hdtPJdiscr,hdtLJnsvtx,hdtPJnsvtx,hdtnumtag,hdtLJsvtxdls,hdtPJsvtxdls,hdtLJsvtxm,hdtPJsvtxm,hdtLJsvtxm_ntrkGT2,hdtPJsvtxm_ntrkGT2});
  SetB({hdtLJdiscr_ntrkGT2,hdtPJdiscr_ntrkGT2,hdtLJdiscr,hdtPJdiscr,hdtLJnsvtx,hdtPJnsvtx,hdtnumtag,hdtLJsvtxdls,hdtPJsvtxdls,hdtLJsvtxm,hdtPJsvtxm,hdtLJsvtxm_ntrkGT2,hdtPJsvtxm_ntrkGT2});
  
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
  DrawCompare(hdtPJdiscr,hsmcPJdiscr);
  plotymin = 9999;
  plotymax  = pp? 1E8:1E6;
  DrawCompare(hdtPJnsvtx,hsmcPJnsvtx);
  plotymax  =9999;
  DrawCompare(hdtPJsvtxdls,hsmcPJsvtxdls);
  // plotylog = false;
  DrawCompare(hdtPJsvtxm,hsmcPJsvtxm);

  // plotylog = true;
  // plotthirdline += "; # s.v. tracks>2";
  // DrawCompare(hdtPJdiscr_ntrkGT2,hsmcPJdiscr_ntrkGT2);
  // plotylog = false;
  // DrawCompare(hdtPJsvtxm_ntrkGT2,hsmcPJsvtxm_ntrkGT2);

}


void taggingplots()
{
  macro m("taggingplots",true);
  
  float fr = 1.0;

  makeplots(-1,-1,fr);
  makeplots(0,20,fr);
  makeplots(20,60,fr);
  makeplots(60,200,fr);

}
