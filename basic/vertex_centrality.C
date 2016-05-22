#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"


void makeplots(int binMin, int binMax)
{
  bool pp = binMin==-1 && binMax==-1; //pp = false, PbPb = true

  auto fmcqcd = config.getfile_djt(pp ? "mcppqcd":"mcPbqcd");
  // auto fmcbfa = config.getfile_djt(pp ? "mcppbfa":"mcPbbfa");
  auto fdtbjt = config.getfile_djt(pp ? "dtppjpf":"dtPbj60");

  auto fdtBJT = config.getfile_djt(pp ? "dtppjpf":"dtPbbjt");


  buildNamesuffix = pp? "pp" : Form("_bin%d%d",binMin, binMax);
  seth(30,-15,15);
  auto hmcvtx = geth("hmcvtx","MC;vz [cm]");
  auto hdtvtx = geth("hdtvtx","Data;vz [cm]");

  seth(50,0,200);
  auto hmcbin = geth("hmcbin","MC;centrality bin");
  auto hdtbin = geth("hdtbin","Data;centrality bin");
  auto hdtbinbjt = geth("hdtbinbjt","Data b-jets;centrality bin");



  // Fill(fmcbfa,{"weight","pthat","bProdCode","vz","bin","jtpt1","refpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL","subidSL","refptSL","pairCodeSL1",
  //                 "refparton_flavorForBSL","refparton_flavorForB1","numTagged",
  //                  "nsvtx1","svtxdls1","svtxm1","nsvtxSL","svtxdlsSL","svtxmSL","svtxntrk1"},[&] (dict &d) {

  //   if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForBSL"])!=5) return; //use bjt sample with b on either side 

  //   fillmchistsfunc(d);

  // });

  Fill(fmcqcd,{"weight","pthat","vz","bin","jtpt1","jtpt2","dphi21"},[&] (dict &d) {
    // if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;
    if (d["pthat"]<pthatcut) return;

    float w = d["weight"];//pp ? d["weight"] : weight1SLPbPb(d);

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>2.) {
      hmcvtx->Fill(d["vz"],w);
      hmcbin->Fill(d["bin"],w);
    }

  });

  Fill(fdtbjt,{"weight","vz","bin","jtpt1","jtpt2","dphi21"},[&] (dict &d) {
    // if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];//pp ? d["weight"] : weight1SLPbPb(d);

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>2.) {
      hdtvtx->Fill(d["vz"],w);
      hdtbin->Fill(d["bin"],w);
    }

    });

  Fill(fdtBJT,{"weight","vz","bin","jtpt1","jtpt2","dphi21","jtptSL","dphiSL1","discr_csvV1_1"},[&] (dict &d) {
    // if (!pp && (d["bin"]<binMin || d["bin"]>=binMax)) return;

    float w = d["weight"];//pp ? d["weight"] : weight1SLPbPb(d);

    if (d["jtpt1"]>100 && d["discr_csvV1_1"]>0.9 && d["jtptSL"]>40 && d["dphiSL1"]>2.) {
      hdtbinbjt->Fill(d["bin"],w);
    }

    });



  SetInc({hdtvtx,hdtbin,hmcvtx,hmcbin});
  SetB({hdtbinbjt});
  SetMC({hdtbinbjt});
  SetData({hdtvtx,hdtbin});
  SetHist({hmcvtx,hmcbin});
  hmcvtx->SetLineWidth(1);
  hmcbin->SetLineWidth(1);

  
  
  NormalizeAllHists();

  // plotylog = true;
  // plotymax = 0.25;
  plotytitle = "Event fractions";
  // aktstring = "anti-k_{T} Pu R=0.4 |#eta|<2.0";  
  aktstring = pp ? "anti-k_{T} R=0.4 PF" : "anti-k_{T} Pu R=0.4 PF";// : Form("PbPb %d-%d %%",binMin/2,binMax/2);
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";

// plotlegenddx = -0.05;
  DrawCompare(hdtvtx,hmcvtx);
  if (!pp) DrawCompare(hdtbin,hmcbin);//,"",hdtbinbjt);

}


void vertex_centrality()
{
  macro m("vertex_centrality",true);
  
  makeplots(-1,-1);
  makeplots(0,200);
  // makeplots(20,60);
  // makeplots(60,200);

}
