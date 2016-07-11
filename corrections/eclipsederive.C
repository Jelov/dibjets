#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"
#include "TProfile.h"

vector<int> binbounds = {0,5,10,15,20,30,40,50,60,200};//{0,20,60,200};//{0,5,10,15,20,30,40,50,60,200};
vector<float> binmean;
vector<TProfile *> profs;
vector<TF1 *>fs;

vector<float> NSfrac;

int getbin(int bin)
{
  int b = 0;
  while (bin>=binbounds[b] && b<(int)binbounds.size()) b++;
  return b-1;
}

void deriveNSfractions(bool data, int mode, bool bjet = false)
{
  if (mode==0) {
    for (unsigned i=0;i<binbounds.size()-1;i++)
      NSfrac.push_back(1);
    return;
  }
  if (mode==2) {
    for (unsigned i=0;i<binbounds.size()-1;i++)
      NSfrac.push_back(0.7);
    return;
  }
  if (mode==4) {
    for (unsigned i=0;i<binbounds.size()-1;i++)
      NSfrac.push_back(1.3);
    return;
  }

  TString fnpp;
  if (bjet)
    fnpp = data ? "dtppjpf" : "mcppbfa";
  else
    fnpp = data ? "dtppjpf" : "mcppqcd";
  auto filepp = config.getfile_djt(fnpp);

  seth(1,0,PI);
  auto hppNS = geth("hppNS");
  auto hppAS = geth("hppAS");

  Fill(filepp,[&] (dict &d) {
    if (bjet && !(d["discr_csvV1_1"]>0.9 && d["discr_csvV1_2"]>0.9)) return;
    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]<PI13) hppNS->Fill(d["dphi21"],d["weight"]);
    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) hppAS->Fill(d["dphi21"],d["weight"]);
  });

  float a_pp = hppNS->Integral()/hppAS->Integral();
  cout<<"a_pp "<<a_pp<<endl;


  TString fnPbPb;
  if (bjet)
    fnPbPb = data ? "dtPbbjt" : "mcPbbfa";
  else
    fnPbPb = data ? "dtPbjcl" : "mcPbqcd";
  auto file = config.getfile_djt(fnPbPb);
  vector<TString> histn;
  for (unsigned i=0;i<binbounds.size()-1;i++)
    histn.push_back(Form("PbPb%d%d",binbounds[i],binbounds[i+1]));
  setv(histn);

  auto hPbPbNS = getv("hPbPbNS");
  auto hPbPbAS = getv("hPbPbAS");

  auto hPbPbNSsig = getv("hPbPbNSsig");

  Fill(file,[&] (dict &d) {
    int b = getbin(d["bin"]);
    float w = d["weight"];
    if (bjet && !(d["discr_csvV1_1"]>0.9 && d["discr_csvV1_2"]>0.9)) return;  
    if (bjet && d["pairCode21"]==0) w*=processweight((int)d["bProdCode"]);
    if (!data && d["pthat"]<50) return;  
    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]<PI13) hPbPbNS[b]->Fill(d["dphi21"],d["weight"]);
    if (d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]>PI23) hPbPbAS[b]->Fill(d["dphi21"],d["weight"]);

    if (!data && d["jtpt1"]>pt1cut && d["jtpt2"]>pt2cut && d["dphi21"]<PI13 && (d["subid2"]==0 && d["refpt2"]>20))
      hPbPbNSsig[b]->Fill(d["dphi21"],d["weight"]);

  });

  vector<double> binmean;
  vector<double> nsfractrue;
  vector<double> nsfracest;

  for (unsigned i=0;i<binbounds.size()-1;i++) {
    float a_PbPb = hPbPbNS[i]->Integral()/hPbPbAS[i]->Integral();
    float x = (1-a_pp/a_PbPb)/(1-a_pp);
    NSfrac.push_back(mode==1 ? x : -1);// = old stupid adding of the signal : 2-x);
    // cout<<binbounds[i]<<" - "<<binbounds[i+1]<<" = "<<NSfrac[i]<<endl;

    nsfractrue.push_back(1-hPbPbNSsig[i]->Integral()/hPbPbNS[i]->Integral());
    nsfracest.push_back(NSfrac[NSfrac.size()-1]);
    binmean.push_back((binbounds[i]+binbounds[i+1])/2);
  }

  auto c = getc();
  TGraph *g = new TGraph(binbounds.size()-1,&binmean[0],&nsfractrue[0]);
  g->Draw();
  SavePlot(c,bjet ? "bjetfractionofNS" : "incfractionofNS");

  auto c2 = getc();
  TGraph *g2 = new TGraph(binbounds.size()-1,&binmean[0],&nsfracest[0]);
  g2->Draw();
  SavePlot(c2,bjet ? "bjetfractionofNSdata" : "incfractionofNSdata");



///////////////////////////// DATA:
// bin<20 : 0.93
// bin>=20 && bin<60 : 0.84
// bin>=60 : 0.39
////////////////////////////


}

void derivefromprofile()
{

  auto file = config.getfile_djt("mcPbqcd");
  auto nt = (TTree *)file->Get("nt");

  for (unsigned i=1;i<binbounds.size();i++) {
    int b1 = binbounds[i-1];
    int b2 = binbounds[i];
    auto p = new TProfile(Form("p%d%d",b1,b2),Form("prof"),100,0,200);

    nt->Project(p->GetName(),"(subid2 == 0 && refpt2 > 20):jtptSignal2",Form("weight*(jtpt1>100&&bin>=%d && bin<%d)",b1,b2));
    profs.push_back(p);

    auto f = new TF1(Form("f%d%d",b1,b2),"exp(-[0]*exp(-[1]*x))");
    f->SetParameters(100,0.1);
    p->Fit(f);
    fs.push_back(f);

    float median = -1/f->GetParameter(1)*log(-1/f->GetParameter(0)*log(0.5));

    auto c = getc();
      TLatex *Tl = new TLatex();
    p->SetMinimum(0);p->SetMaximum(1);
    p->Draw();
    f->Draw("same");
    Tl->DrawLatexNDC(0.6,0.4,Form("PbPb bin %d-%d",b1,b2));
    Tl->DrawLatexNDC(0.6,0.35,Form("median = %.2f",median));

    TLine *l1 = new TLine(median,0,median, f->Eval(median));
    l1->Draw();
    TLine *l2 = new TLine(0,0.5,median, f->Eval(median));
    l2->Draw();

    SavePlot(c,Form("fit%d%d",b1,b2));
  }

}

float pt = 40;
vector<float> prob;
vector<float> meanb;

void ScaleVisibleBins(TH1F *h, float SF)
{
  float underv = h->GetBinContent(0);
  float undere = h->GetBinError(0);

  float i = h->Integral();
  h->Scale(SF);

  h->SetBinContent(0,underv+i*(1-SF));
  h->SetBinError(0,undere+i*(1-SF));
}


void derivefromNS(bool data = false)
{
  // float shift = mode == 2 ? -2 : mode*2; //mode = 0,1,2 shift = 0,2,-2
  // cout<<"shift = "<<shift<<endl;

  auto file = config.getfile_djt(data ? "dtPbjcl" : "mcPbqcd");
  auto nt = (TTree *)file->Get("nt");

  for (unsigned i=1;i<binbounds.size();i++) {
    int b1 = binbounds[i-1];
    int b2 = binbounds[i];
    seth(71,38,180); //71,38
    auto h = geth(Form("h%d%d",b1,b2)); //allowing one more bin for overflow
    seth(b2-b1,b1,b2);
    auto hb = geth(Form("hb%d%d",b1,b2));

   
    TString mcappendix = data ? "" : "&& pthat>50";//"&& subid2!=0 && pthat>50"; //"&& pthat>50";//"&& !(subid2==0 && refpt2>20) && pthat>50"; //"&& pthat>80";//
   
    //Form("jtpt2+%f",shift)
    nt->Project(h->GetName(),"jtpt2", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05 %s)",b1,b2,mcappendix.Data()));
    nt->Project(hb->GetName(),"bin", Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphi21<1.05 %s)",b1,b2,mcappendix.Data()));

    ScaleVisibleBins(h,NSfrac[i-1]);

    h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));

    // auto p = new TProfile(Form("p%d%d",b1,b2),Form("prof"),71,38,180);
    // nt->Project(p->GetName(),"(subid2 == 0 && refpt2 > 20):jtptSignal2",Form("weight*(jtpt1>100&&bin>=%d && bin<%d && dphiSignal21<1.05 && !(subid2==0 && refpt2>20) && pthat>80)",b1,b2));


    auto g = getCDFgraph(h);
    g->GetXaxis()->SetTitle("p_{T,2} threshold [GeV]");
    g->GetYaxis()->SetTitle("found fraction");

    auto gtemp = getCDFgraph(h);
    meanb.push_back(hb->GetMean());
    prob.push_back(gtemp->Eval(pt));

    float c0=g->Eval(50);
    cout<<"prob at 50: "<<c0<<endl;

    auto gf = new TFile("graph.root","recreate");
    gtemp->Write();
    gf->Close();

    // auto f = new TF1(Form("f%d%d",b1,b2),"1-[0]*exp(-[1]*(x-40))",40,180);
    // f->SetParameters(0.1,0.1);
    // // f->FixParameter(1,0.08);

    // auto f = new TF1(Form("f%d%d",b1,b2),"TMath::Erf((x-[0])/[1])",40,180);
    // f->SetParameters(40,10);
    // f->FixParameter(1,25);

    auto f = new TF1(Form("f%d%d",b1,b2),"exp(-[0]*exp(-[1]*x))",40,180);
    f->SetParameters(100,0.1);
    // f->FixParameter(1,0.11); //!!!!!!!!!!

    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    g->Fit(f,"RM");
    fs.push_back(f);
    binmean.push_back(hb->GetMean());

    float median = -1/f->GetParameter(1)*log(-1/f->GetParameter(0)*log(0.5));

    Draw({h});

    // h->Rebin(2);
    auto gcoarse = getCDFgraph(h);

    auto c = getc();
          TLatex *Tl = new TLatex();
    gcoarse->SetMinimum(0);g->SetMaximum(1);
    gcoarse->Draw("AP");
    f->Draw("same");
    Tl->DrawLatexNDC(0.6,0.55,"y=e^{-a e^{-b x} }");
    Tl->DrawLatexNDC(0.6,0.50,Form("a = %.1f",f->GetParameter(0)));
    Tl->DrawLatexNDC(0.6,0.45,Form("b = %.2f",f->GetParameter(1)));
    Tl->DrawLatexNDC(0.6,0.4,Form("PbPb bin %d-%d",b1,b2));
    Tl->DrawLatexNDC(0.6,0.35,Form("median = %.2f",median));

    TLine *l1 = new TLine(median,0,median, f->Eval(median));
    l1->Draw();
    TLine *l2 = new TLine(0,0.5,median, f->Eval(median));
    l2->Draw();

    // p->SetMarkerColor(kRed);
    // p->SetMarkerSize(1);
    // p->Draw("same");
    SavePlot(c,Form("fit%d%d",b1,b2));//return;
  }

}

void eclipsederive(bool data = false, int mode = 0)
{
  TString mcdt(data ? "dt" : "mc");

  macro m(Form("eclipsederive0707_%s_%d",mcdt.Data(),mode));

  // //just to know NS fractions in b-jet data
  // deriveNSfractions(data,mode,true);
  // cout<<"NSfracinc"<<mcdt<<" = {";
  // for (unsigned i=0;i<binbounds.size()-1;i++) {
  //   cout<<NSfrac[i];
  //   if (i<binbounds.size()-2) cout<<",";
  // }
  // cout<<"};"<<endl;

  deriveNSfractions(data,mode);

  derivefromNS(data);

  //derivefromprofile();

  auto ggg = new TGraph(meanb.size(),&meanb[0],&prob[0]);
  auto cc = getc();
  ggg->Draw("AP");
  SavePlot(cc,"probvsbin");



  vector<float> x = binmean;
  vector<float> xerr;
  vector<float> y,yerr;
  vector<float> z,zerr;


  for (unsigned i=1;i<binbounds.size();i++) {
    // int b1 = binbounds[i-1];
    // int b2 = binbounds[i];
    // x.push_back((b1+b2)/2);

    xerr.push_back(0);
    y.push_back(fs[i-1]->GetParameter(0));
    yerr.push_back(fs[i-1]->GetParError(0));

    z.push_back(fs[i-1]->GetParameter(1));
    zerr.push_back(fs[i-1]->GetParError(1));
  }
  auto c2 = getc();
  auto g= new TGraphErrors(x.size(),&x[0],&y[0],&xerr[0],&yerr[0]);
  g->Draw();
  SavePlot(c2,"par0cent");

  auto c3 = getc();
  auto g2 = new TGraphErrors(x.size(),&x[0],&z[0],&xerr[0],&zerr[0]);
  g2->Draw();
  SavePlot(c3,"par1cent");


  std::ofstream res (plotfoldername+"/eclipserapams.txt", std::ofstream::out);

  res<<"binmean"<<mcdt<<" = {";
  for (unsigned i=0;i<binmean.size();i++) {
    res<<binmean[i];
    if (i<binmean.size()-1) res<<",";
  }
  res<<"};"<<endl;

  res<<"NSfracinc"<<mcdt<<" = {";
  for (unsigned i=0;i<binbounds.size()-1;i++) {
    res<<NSfrac[i];
    if (i<binbounds.size()-2) res<<",";
  }
  res<<"};"<<endl;

  res<<"par0"<<mcdt<<" = {";
  for (unsigned i=1;i<binbounds.size();i++) {
    int b1 = binbounds[i-1];
    int b2 = binbounds[i];
    res<<fs[i-1]->GetParameter(0);
    if (i<binbounds.size()-1) res<<",";
  }
  res<<"};"<<endl;

  res<<"par1"<<mcdt<<" = {";
  for (unsigned i=1;i<binbounds.size();i++) {
    int b1 = binbounds[i-1];
    int b2 = binbounds[i];
    res<<fs[i-1]->GetParameter(1);
    if (i<binbounds.size()-1) res<<",";
  }
  res<<"};"<<endl;

  res.close();

}