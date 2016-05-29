#include "../helpers/plotting.h"
#include "../helpers/config.h"
#include "TProfile.h"

void drawnearside()
{
  macro m("drawnearside");

  auto fdt = config.getfile_djt("dtPbj60");
  auto ntdt = (TTree *)fdt->Get("nt");

  auto fmc = config.getfile_djt("mcPbqcd");
  auto ntmc = (TTree *)fmc->Get("nt");

  int b1=0;
  int b2=5;

  seth(35,40,180);
  auto hdt = geth(Form("hdt%d%d",b1,b2),"Data; Max jet p_{T} [GeV]; Event fractions");
  auto hmc = geth(Form("hmc%d%d",b1,b2),"MC; Max jet p_{T} [GeV]; Event fractions");
 
  ntdt->Project(hdt->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05)",b1,b2));
  ntmc->Project(hmc->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05 && pthat>50)",b1,b2));


  auto gdt = getCDF(hdt);
  auto gmc = getCDF(hmc);

  SetMC({hmc});
  SetData({hdt});

  Normalize({hdt,hmc});
  aktstring = "#Delta#phi<1/3#pi";
  plotsecondline = "PbPb "+FloatToStr(b1/2.)+"-"+FloatToStr(b2/2.)+"%";
  plottextposy = 0.5;
  plotymax = 0.2;
  Draw({hdt,hmc});




  gdt->SetLineColor(darkred);
  gmc->SetLineColor(darkgreen);
  gdt->SetLineWidth(3);
  gmc->SetLineWidth(3);

  gdt->SetMinimum(0);
  gdt->SetMaximum(1);

  gdt->GetXaxis()->SetTitle("p_{T,2} threshold [GeV]");
  gdt->GetYaxis()->SetTitle("found fraction"); //c.d.f. of near-side jets

  auto c = getc();
  gdt->Draw("hist");
  gmc->Draw("hist,same");
  plotlegendpos = BottomRight;
  // plotlegenddx = -0.1;
  auto l=getLegend();
  l->AddEntry(gdt,"Data","L");
  l->AddEntry(gmc,"MC","L");

  l->SetHeader("PbPb "+FloatToStr(b1/2.)+"-"+FloatToStr(b2/2.)+"%");
  l->Draw();
  SavePlots(c,"cdfDataMC");

}