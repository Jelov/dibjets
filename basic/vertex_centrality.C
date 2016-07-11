#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "../corrections/eclipseclosure.C"

void makeplots(int binMin, int binMax)
{
  bool pp = binMin==-1 && binMax==-1; //pp = false, PbPb = true

  auto fmcqcd = config.getfile_djt(pp ? "mcppqcd":"mcPbqcd");
  auto fmcbfa = config.getfile_djt(pp ? "mcppbfa":"mcPbbfa");
  auto fdtinc = config.getfile_djt(pp ? "dtppjpf":"dtPbj60");
  auto fdtBJT = config.getfile_djt(pp ? "dtppjpf":"dtPbbjt");


  buildNamesuffix = pp? "pp" : Form("_bin%d%d",binMin, binMax);
  seth(30,-15,15);
  auto hmcvtx = geth("hmcvtx","MC;vz [cm]");
  auto hdtvtx = geth("hdtvtx","Data;vz [cm]");

  seth(50,0,200);
  auto hmcbin = geth("hmcbin","MC Inc;centrality bin");
  auto hdtbin = geth("hdtbin","Data Inc;centrality bin");
  auto hdtbinbjt = geth("hdtbinbjt","Data b-jets AS raw;centrality bin");


  auto hmcbinASraw = geth("hmcbinASraw","MC Inc AS raw;centrality bin");
  auto hdtbinASraw = geth("hdtbinASraw","Data Inc AS raw;centrality bin");

  auto hmcbinAS = geth("hmcbinAS","MC Inc Subtracted;centrality bin");
  auto hdtbinAS = geth("hdtbinAS","Data Inc Subtracted;centrality bin");
  auto hdtbinbjtAS = geth("hdtbinbjtAS","Data b-jets Subtracted;centrality bin");

  auto hmcbinNS = geth("hmcbinNS","MC Inc NS;centrality bin");
  auto hdtbinNS = geth("hdtbinNS","Data Inc NS;centrality bin");
  auto hdtbinbjtNS = geth("hdtbinbjtNS","Data b-jets NS;centrality bin");


  auto hmcbjtbinASraw = geth("hmcbjtbinASraw","MC bjet AS raw;centrality bin");
  auto hmcbjtbinAStrue = geth("hmcbjtbinAStrue","MC bjet AS true;centrality bin");
  auto hmcbjtbinAS = geth("hmcbjtbinAS","MC bjet AS Subtracted;centrality bin");
  auto hmcbjtbinNS = geth("hmcbjtbinNS","MC bjet NS;centrality bin");


  Fill(fmcbfa,{"weight","pthat","vz","bin","jtpt1","jtptSL","dphiSL1","pairCodeSL1","bProdCode","subidSL","refptSL","discr_csvV1_1",
                "refparton_flavorForB1","refparton_flavorForBSL","numTagged"},[&] (dict &d) {
    if (d["numTagged"]>6) return;
    if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForBSL"])!=5) return; 

    if (d["pthat"]<pthatcut) return;
    float w = d["weight"];//pp ? weight1SLpp(d) : weight1SLPbPb(d);

    if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]>PI23) {
      hmcbjtbinASraw->Fill(d["bin"],w);
      hmcbjtbinAS->Fill(d["bin"],w);

      if (d["pairCodeSL1"]==0)
        hmcbjtbinAStrue->Fill(d["bin"],w);
    }

    if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]<PI13)
      hmcbjtbinNS->Fill(d["bin"],w);


  });


  //   Fill(fmcqcd,{"weight","pthat","vz","bin","jtpt1","jtptSL","dphiSL1","pairCodeSL1","bProdCode","subidSL","refptSL","discr_csvV1_1",
  //               "refparton_flavorForB1","refparton_flavorForBSL","numTagged"},[&] (dict &d) {
  //   if (d["numTagged"]>6) return;
  //   if (abs(d["refparton_flavorForB1"])==5 || abs(d["refparton_flavorForBSL"])==5) return; 

  //   if (d["pthat"]<pthatcut) return;
  //   float w = d["weight"];//pp ? weight1SLpp(d) : weight1SLPbPb(d);

  //   if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]>PI23) {
  //     hmcbjtbinASraw->Fill(d["bin"],w);
  //     hmcbjtbinAS->Fill(d["bin"],w);

  //     if (d["pairCodeSL1"]==0)
  //       hmcbjtbinAStrue->Fill(d["bin"],w);
  //   }

  //   if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]<PI13)
  //     hmcbjtbinNS->Fill(d["bin"],w);


  // });

  // Fill(fmcbfa,{"weight","pthat","vz","bin","jtpt1","jtpt2","dphi21","pairCode21","bProdCode","subid2","refpt2","discr_csvV1_1",
  //               "refparton_flavorForB1","refparton_flavorForB2","numTagged"},[&] (dict &d) {
  //   if (d["numTagged"]>6) return;
  //   if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForB2"])!=5) return; 

  //   if (d["pthat"]<pthatcut) return;
  //   float w = d["weight"];//pp ? weight1SLpp(d) : weight1SLPbPb(d);

  //   if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>PI23) {
  //     hmcbjtbinASraw->Fill(d["bin"],w);

  //     if (d["pairCode21"]==0)
  //       hmcbjtbinAStrue->Fill(d["bin"],w);
  //   }
  // });







  // Fill(fmcqcd,{"weight","pthat","vz","bin","jtpt1","jtptSL","dphiSL1","pairCodeSL1","bProdCode","subidSL","refptSL","discr_csvV1_1",
  //               "refparton_flavorForB1","refparton_flavorForBSL","numTagged"},[&] (dict &d) {
  //   if (d["numTagged"]>6) return;
  //   if (abs(d["refparton_flavorForB1"])==5 || abs(d["refparton_flavorForBSL"])==5) return; 

  //   if (d["pthat"]<pthatcut) return;
  //   float w = pp ? weight1SLpp(d) : weight1SLPbPb(d);

  //   if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]>PI23) {
  //     hmcbjtbinASraw->Fill(d["bin"],w);

  //   if (d["pairCodeSL1"]==0)
  //     hmcbjtbinAStrue->Fill(d["bin"],w);
  //   }
  // });




  Fill(fmcqcd,{"weight","pthat","vz","bin","jtpt1","jtpt2","dphi21"},[&] (dict &d) {
    // if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;
    if (d["pthat"]<pthatcut) return;

    float w = d["weight"];//pp ? d["weight"] : weight1SLPbPb(d);
    float ew = eclipseWeightmc(d["jtpt2"],d["bin"]);

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>PI23) {
      hmcvtx->Fill(d["vz"],w);
      hmcbin->Fill(d["bin"],w);

      hmcbinASraw->Fill(d["bin"],w);
      hmcbinAS->Fill(d["bin"],w*ew);

    }

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]<PI13) {
      hmcbinNS->Fill(d["bin"],w*ew);

    }
  });

  Fill(fdtinc,{"weight","vz","bin","jtpt1","jtpt2","dphi21"},[&] (dict &d) {
    float w = d["weight"];
    float ew = eclipseWeightdt(d["jtpt2"],d["bin"]);

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>PI23) {
      hdtvtx->Fill(d["vz"],w);
      hdtbin->Fill(d["bin"],w);

      hdtbinASraw->Fill(d["bin"],w);
      hdtbinAS->Fill(d["bin"],w*ew);
    }

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]<PI13) {
      hdtbinNS->Fill(d["bin"],w*ew);

    }

    });

  Fill(fdtBJT,{"weight","vz","bin","jtpt1","jtpt2","dphi21","jtptSL","dphiSL1","discr_csvV1_1","discr_csvV1_2","numTagged"},[&] (dict &d) {
    if (d["numTagged"]>6) return;
    float w = d["weight"];


    if (d["jtpt1"]>100 && d["jtptSL"]>40 && d["dphiSL1"]>PI23 && d["discr_csvV1_1"]>0.9) {
      hdtbinbjt->Fill(d["bin"],w);
      hdtbinbjtAS->Fill(d["bin"],w);
    }

    if (d["jtpt1"]>100 && d["jtptSL"]>40 && d["dphiSL1"]<PI13 && d["discr_csvV1_1"]>0.9) {
      hdtbinbjtNS->Fill(d["bin"],w);
    }

    // if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>PI23 && d["discr_csvV1_1"]>0.9 && d["discr_csvV1_2"]>0.9) {
    //   hdtbinbjt->Fill(d["bin"],w);
    //   hdtbinbjtAS->Fill(d["bin"],w);
    // }

    // if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]<PI13 && d["discr_csvV1_1"]>0.9 && d["discr_csvV1_2"]>0.9) {
    //   hdtbinbjtNS->Fill(d["bin"],w);
    // }
    
    });

  hmcbinAS->Add(hmcbinNS,-1);
  hdtbinAS->Add(hdtbinNS,-1);
  hdtbinbjtAS->Add(hdtbinbjtNS,-1);
  hmcbjtbinAS->Add(hmcbjtbinNS,-1);


 
 

  SetInc({hdtvtx,hdtbin,hmcvtx,hmcbin,hmcbinAS,hdtbinAS});
  SetB({hdtbinbjt,hmcbjtbinASraw,hmcbjtbinAStrue,hdtbinbjtAS,hmcbjtbinAS});
  SetMC({hdtbinbjt,hmcbinAS,hmcbjtbinASraw});
  SetData({hdtvtx,hdtbin,hdtbinAS,hmcbjtbinAStrue,hdtbinbjtAS,hmcbjtbinAS});
  SetHist({hmcvtx,hmcbin});
  hmcvtx->SetLineWidth(1);
  hmcbin->SetLineWidth(1);

  

  NormalizeAllHists();

  plotymax = 0.2;
  plotratiomin = 0.5;
  plotratiomax = 1.5;

  // if (!pp) Draw({hmcbinAS,hdtbinAS,hdtbinASraw});

  plotdivide = true;
  // plotdiffmax = 
  // plotylog = true;
  // plotymax = 0.25;
  plotytitle = "Event fractions";
  //aktstring = pp ? "anti-k_{T} R=0.4 PF" : "anti-k_{T} Pu R=0.4 PF";// : Form("PbPb %d-%d %%",binMin/2,binMax/2);
  // plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  // plotthirdline = "#Delta#phi>2/3#pi";

// plotlegenddx = -0.05;

  // DrawCompare(hdtvtx,hmcvtx);
  if (!pp) {
    // DrawCompare(hdtbin,hmcbin);//,"",hdtbinbjt);
    DrawCompare(hdtbinAS,hdtbinASraw); //,"",hdtbinbjtAS
    DrawCompare(hmcbinAS,hmcbinASraw);

    DrawCompare(hmcbjtbinASraw,hmcbjtbinAStrue);
    
    DrawCompare(hmcbjtbinASraw,hmcbinASraw);    
    DrawCompare(hdtbinbjt,hdtbinbjtAS);

    DrawCompare(hmcbjtbinASraw,hmcbjtbinAS);

    DrawCompare(hmcbinAS,hmcbjtbinAS);
  }

}


void vertex_centrality()
{
  macro m("vertex_centrality",true);
  
  // makeplots(-1,-1);
  makeplots(0,200);
  // makeplots(20,60);
  // makeplots(60,200);

}
