#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

void drawdphideta()
{
  macro m("drawdphideta_0p95",true);

  csvcut = 0.95;

  seth(25,0,3.142,25,0,4);

  auto hdphidetasig = geth2d("hdphidetasig",";#Delta#phi;#Delta#eta");
  auto hdphidetabkg = geth2d("hdphidetabkg",";#Delta#phi;#Delta#eta");

  auto hdphidetasig12 = geth2d("hdphidetasig12",";#Delta#phi;#Delta#eta");
  auto hdphidetabkg12 = geth2d("hdphidetabkg12",";#Delta#phi;#Delta#eta");


  auto fmcPb = config.getfile_djt("mcPbbfa");

  Fill(fmcPb,[&] (dict &d) {
      if (d["pthat"]<pthatcut) return;
      float w = weight1SLPbPb(d);
      if (d["bin"]>20) return;

      float w2 = d["weight"];
      if (d["pairCodeSignal21"]==0) w2*=processweight((int)d["bProdCode"]);


      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) {
        if (IsSignal(d))
          hdphidetasig->Fill(d["dphiSL1"],abs(d["jteta1"]-d["jtetaSL"]),w);
        else
          hdphidetabkg->Fill(d["dphiSL1"],abs(d["jteta1"]-d["jtetaSL"]),w);
      }

      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtptSignal2"]>pt2cut && d["discr_csvV1_Signal2"]>csvcut)
          hdphidetasig12->Fill(d["dphiSignal21"],abs(d["jteta1"]-d["jtetaSignal2"]),w2);

      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtpt2"]>pt2cut && d["discr_csvV1_2"]>csvcut && !IsSignal(d["subid2"],d["refpt2"]))
          hdphidetabkg12->Fill(d["dphi21"],abs(d["jteta1"]-d["jteta2"]),w);


      // if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtpt2"]>pt2cut) {
      //   if (IsSignal(d["subid2"],d["refpt2"]))
      //     hdphidetasig12->Fill(d["dphi21"],abs(d["jteta1"]-d["jteta2"]),w2);
      //   else
      //     hdphidetabkg12->Fill(d["dphi21"],abs(d["jteta1"]-d["jteta2"]),w);
      // }

    
      });

  auto c1 = getc();
  c1->SetLogz();
  hdphidetasig->SetMaximum(5E-9);
  hdphidetasig->Draw("colz");
  SavePlots(c1,"signaldphideta");

  auto c2 = getc();
  c2->SetLogz();
  hdphidetabkg->SetMaximum(5E-9);
  hdphidetabkg->Draw("colz");
  SavePlots(c2,"bkgdphideta");


  auto c3 = getc();
  c3->SetLogz();
  hdphidetasig12->SetMaximum(5E-9);//1E-7);
  hdphidetasig12->Draw("colz");
  SavePlots(c3,"signaldphideta12");

  auto c4 = getc();
  c4->SetLogz();
  hdphidetabkg12->SetMaximum(5E-9);//1E-7);
  hdphidetabkg12->Draw("colz");
  SavePlots(c4,"bkgdphideta12");

}
