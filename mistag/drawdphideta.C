#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

void drawdphideta()
{
  macro m("drawdphideta");

  buildh(50,0,3.142,50,0,4);

  auto hdphidetasig = geth2d("hdphidetasig");
  auto hdphidetabkg = geth2d("hdphidetabkg");


  auto fmcPb = config.getfile_djt("mcPbbfa");

  Fill(fmcPb,{"pthat","weight","jtpt1","refpt1","bProdCode","jtptSL","refptSL","dphiSL1","refparton_flavorForB1","subidSL","bin","pairCodeSL1","jteta1","jtetaSL"},[&] (dict d) {
      if (d["pthat"]<pthatcut) return;
      float w = weight1SLPbPb(d);

      if (d["jtpt1"]>pt1cut && d["refpt1"]>50 && abs(d["refparton_flavorForB1"])==5 && d["jtptSL"]>pt2cut) {
        if (IsSignal(d))
          hdphidetasig->Fill(d["dphiSL1"],abs(d["jteta1"]-d["jtetaSL"]),w);
        else
          hdphidetabkg->Fill(d["dphiSL1"],abs(d["jteta1"]-d["jtetaSL"]),w);
      }





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

}
