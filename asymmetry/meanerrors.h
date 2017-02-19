#include "../helpers/plotting.h"

float integralerror(TH1F *h)
{
  double e;
  h->IntegralAndError(1,h->GetNbinsX(),e);
  return e;
}

void getsubmeanerror(TH1F *h1, TH1F *h2, float &m, float &e)
{
  float m1 = h1->GetMean();
  float e1 = h1->GetMeanError();
  float s1 = h1->GetStdDev();
  double i1 = h1->Integral();
  float N1 = h1->Integral();
  double ie1 = integralerror(h1);
  
  float m2 = h2->GetMean();
  float e2 = h2->GetMeanError();
  float s2 = h2->GetStdDev();
  float i2 = h2->Integral();
  float N2 = h2->Integral();
  double ie2 = integralerror(h2);  

  float a1 = i1/(i1-i2);
  float a2 = i2/(i1-i2);
   
  m = m1*a1-m2*a2;
  
  float dmdN1 = N2*(m2-m1)/((N1-N2)*(N1-N2));
  float dmdN2 = N1*(m1-m2)/((N1-N2)*(N1-N2));
  float dmdm1 = N1/(N1-N2);
  float dmdm2 = -N2/(N1-N2);
  
  float err2 = N1*dmdN1*dmdN1+N2*dmdN2*dmdN2+e1*e1*dmdm1*dmdm1+e2*e2*dmdm2*dmdm2;
  // float err2 = ie1*ie1*dmdN1*dmdN1+ie2*ie2*dmdN2*dmdN2+e1*e1*dmdm1*dmdm1+e2*e2*dmdm2*dmdm2;
  e = sqrt(err2);
}

void generateunif(TH1F *h, TH1F *h2, int nev1, int nev2)
{
  TF1 f1("f1","pol0",0,5);
  TF1 f2("f2","pol0",5,10);
  f1.SetParameter(0,1);
  f2.SetParameter(0,2);
  
  // for gauss
  // f1.SetParameters(100,5,3);
  // f2.SetParameters(100,4,3);
  


  h->FillRandom("f1",nev1);
  
  for (int i=0;i<nev2;i++) {
    float r = f2.GetRandom();
    if (r>5)
    h2->Fill(r,0.6);
  }
  

}

void generategaus(TH1F *h, TH1F *h1, TH1F *h2, int nev1, int nev2)
{
  TF1 f1("f1","gaus(3)",0,10);
  TF1 f2("f2","gaus(3)",0,10);

  f1.SetParameters(100,10,2);
  f2.SetParameters(100,4,2);
  
  // h->FillRandom("f1",nev1);
  // h1->FillRandom("f1",nev1);
  // h2->FillRandom("f2",nev2);  

  for (int i=0;i<nev1;i++) {
    float r = f1.GetRandom();
    h->Fill(r);
    h1->Fill(r);
  }

  for (int i=0;i<nev2;i++) {
    float r = f2.GetRandom();
    float r2 = f2.GetRandom();
    h->Fill(r,2);
    h2->Fill(r2,2);
  }

  cout<<"True :"<<h1->GetMean()<<" ± "<<h1->GetMeanError()<<endl;

}

void meanerrors()
{
  int N=10000;
  seth(10,7,9);
  auto htruemean = geth("htruemean");
  auto hminemean = geth("hminemean");
  auto hrootmean = geth("hrootmean");

  seth(10,0,1);
  auto htrueerr = geth("htrueerr");
  auto hmineerr = geth("hmineerr");
  auto hrooterr = geth("hrooterr");


  gRandom->SetSeed(0);

for (int i=0;i<N;i++) {
  seth(5,0,10);

  auto h = geth("h");
  auto h1 = geth("h1");
  auto h2 = geth("h2");
  auto hsub = geth("hsub");

  // generateunif(h,h2,100,10);
  generategaus(h, h1,h2,300,300);

  htruemean->Fill(h1->GetMean());
  htrueerr->Fill(h1->GetMeanError());

  float m,e;
  
  getsubmeanerror(h,h2,m,e);
  
  cout<<"Mine :"<<m<<" ± "<<e<<endl;
  
  hminemean->Fill(m);
  hmineerr->Fill(e);



  hsub->Add(h,h2,1,-1);
  
  cout<<"Root :"<<hsub->GetMean()<<" ± "<<hsub->GetMeanError()<<endl;

  hrootmean->Fill(hsub->GetMean());
  hrooterr->Fill(hsub->GetMeanError());


  // h1->SetMinimum(0);
  // h->SetMarkerColor(kRed);
  // h1->SetMarkerColor(kGreen);

  // h->Draw();
  // h->Draw("same");
  // h2->Draw("same");

  hsub->SetMarkerStyle(kOpenSquare);

  plotymin = -50;
  if (N==1) Draw({h,h1,h2,hsub});


  delete h;
  delete h1;
  delete h2;
  delete hsub;

  }

  plotlegendpos = TopLeft;
  plotputmean = true;
  plotusestderror = false;
  Draw({htruemean, hminemean, hrootmean});
  Draw({htrueerr, hmineerr, hrooterr});

}