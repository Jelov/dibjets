#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

// TF1 *f;

void findfunc()
{

  //dijet case
  // auto fmcppbjt = config.getfile_djt("mcppbfa");
  // auto nt = (TTree *)fmcppbjt->Get("nt");
  // seth(18,40,220);
  // auto hd = geth("hd");
  // auto hn = geth("hn");
  // nt->Project("hd","jtpt2","weight*(pthat>80 && jtpt1>100 && abs(refparton_flavorForB1)==5 && abs(refparton_flavorForB2)!=5 && dphi21>2.1)");
  // nt->Project("hn","jtpt2","weight*(pthat>80 && jtpt1>100 && abs(refparton_flavorForB1)==5 && abs(refparton_flavorForB2)!=5 && dphi21>2.1 && discr_csvV1_2>0.9)");
  // hn->Divide(hn,hd,1,1,"B")
  // hn->Fit(f)
  

  auto fpp = new TF1("fpp","expo",40,200);
  auto f1 = new TF1("fPb1","expo",40,200);
  auto f2 = new TF1("fPb2","expo",40,200);
  auto f3 = new TF1("fPb3","expo",40,200);
  seth(18,40,200);


  auto ntpp = (TTree *)config.getfile_inc("mcppqcd")->Get("nt");

  auto hd = geth("hd");
  auto hn = geth("hn");
  ntpp->Project("hd","jtpt","weight*(pthat>50 && refpt>20 && abs(refparton_flavorForB)!=5)");
  ntpp->Project("hn","jtpt","weight*(pthat>50 && refpt>20 && abs(refparton_flavorForB)!=5 && discr_csvV1>0.9)");
  hn->Divide(hn,hd,1,1,"B");
  hn->Fit(fpp);

  auto nt = (TTree *)config.getfile_inc("mcPbqcd")->Get("nt");

  auto hd1 = geth("hd1");
  auto hn1 = geth("hn1");
  nt->Project("hd1","jtpt","weight*(pthat>50 && refpt>20 && bin < 20 && abs(refparton_flavorForB)!=5)");
  nt->Project("hn1","jtpt","weight*(pthat>50 && refpt>20 && bin < 20 && abs(refparton_flavorForB)!=5 && discr_csvV1>0.9)");
  hn1->Divide(hn1,hd1,1,1,"B");
  hn1->Fit(f1);

  auto hd2 = geth("hd2");
  auto hn2 = geth("hn2");
  nt->Project("hd2","jtpt","weight*(pthat>50 && refpt>20 && bin>=20 && bin<60 && abs(refparton_flavorForB)!=5)");
  nt->Project("hn2","jtpt","weight*(pthat>50 && refpt>20 && bin>=20 && bin<60 && abs(refparton_flavorForB)!=5 && discr_csvV1>0.9)");
  hn2->Divide(hn2,hd2,1,1,"B");
  hn2->Fit(f2);

  auto hd3 = geth("hd3");
  auto hn3 = geth("hn3");
  nt->Project("hd3","jtpt","weight*(pthat>50 && refpt>20 && bin>=60 && abs(refparton_flavorForB)!=5)");
  nt->Project("hn3","jtpt","weight*(pthat>50 && refpt>20 && bin>=60 && abs(refparton_flavorForB)!=5 && discr_csvV1>0.9)");
  hn3->Divide(hn3,hd3,1,1,"B");
  hn3->Fit(f3);

  auto fout = new TFile("../ntuplizer/BXmistagfunc.root","recreate");
  fout->cd();
  fpp->Write();
  f1->Write();
  f2->Write();
  f3->Write();
  hn->Write();
  hn1->Write();
  hn2->Write();
  hn3->Write();



  auto c = getc();
  f1->SetLineColor(kRed);
  f2->SetLineColor(kGreen);
  f3->SetLineColor(kOrange);
  fpp->SetLineColor(kBlue);

  auto l = getLegend();
  l->AddEntry(fpp,"pp","L");
  l->AddEntry(f1,"bin<20","L");
  l->AddEntry(f2,"20<bin<60","L");
  l->AddEntry(f3,"bin>60","L");


  f1->Draw();
  f2->Draw("same");
  f3->Draw("same");
  fpp->Draw("same");

  l->Draw();


  c->SaveAs("func.pdf");

  fout->Close();

}


// void run()
// {
//   float pthat = 80;//pthatcut;

//   auto fmcppbjt = config.getfile_djt("mcppbfa");
//   auto fmcppqcd = config.getfile_djt("mcppqcd");
//   auto fdtppbjt = config.getfile_djt("dtppjpf");

//   seth(18,40,220);
//   auto hjt2BXtrue = geth("hjt2BXtrue");
//   auto hjt2BXsim = geth("hjt2BXsim");



//   Fill(fmcppbjt,{}, [&] (dict &d) 
//   {
//     if (d["pthat"]<pthat) return;
//     float w=weight1SLpp(d);
//     if (abs(d["refparton_flavorForB1"])!=5 && abs(d["refparton_flavorForBSL"])!=5) return;

//     if (d["jtpt1"]>pt1cut && abs(d["refparton_flavorForB1"])==5 && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) {
//       if (abs(d["refparton_flavorForB2"])!=5) hjt2BXtrue->Fill(d["jtpt2"],w);
//       float r = gRandom->Uniform();

//     }
 

//   });
// }

void BXmistags()
{
  // macro m("BXmistags",true);
  findfunc();
  // run();
}