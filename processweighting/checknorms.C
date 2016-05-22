#include "../helpers/plotting.h"
float N12_data=0,N13_data=0,N12_fcr =0,N12_fex =0,N12_gsp =0,N13_fcr =0,N13_fex =0,N13_gsp =0;//,N23_data=0,N23_fcr =0,N23_fex =0,N23_gsp=0;

float //N12allmcNS=0, N12fcrNS=0,N12fexNS=0,N12gspNS=0,
      N13allmcNS=0,N13fcrNS=0,N13fexNS=0,N13gspNS=0,
      // N12alldataNS=0,
      N13alldataNS=0,
      N12allmcAS=0,N12fcrAS=0,N12fexAS=0,N12gspAS=0,
      N13allmcAS=0,N13fcrAS=0,N13fexAS=0,N13gspAS=0,
      N12alldataAS=0,N13alldataAS=0;

void loadnumbers2()
{

  auto fNS = new TFile("flavorProcesshists_mcNS.root");
  auto h13allmcNS = (TH1F *)fNS->Get("h13all");
  auto h13fcrNS = (TH1F *)fNS->Get("h13fcr");
  auto h13fexNS = (TH1F *)fNS->Get("h13fex");
  auto h13gspNS = (TH1F *)fNS->Get("h13gsp");
  auto f2NS = new TFile("flavorProcesshists_dataNS.root");
  // auto h12alldataNS = (TH1F *)f2NS->Get("h12all");
  auto h13alldataNS = (TH1F *)f2NS->Get("h13all");


  auto fAS = new TFile("flavorProcesshists_mcAS.root");
  auto h12allmcAS = (TH1F *)fAS->Get("h12all");
  auto h13allmcAS = (TH1F *)fAS->Get("h13all");
  auto h12fcrAS = (TH1F *)fAS->Get("h12fcr");
  auto h12fexAS = (TH1F *)fAS->Get("h12fex");
  auto h12gspAS = (TH1F *)fAS->Get("h12gsp");
  auto h13fcrAS = (TH1F *)fAS->Get("h13fcr");
  auto h13fexAS = (TH1F *)fAS->Get("h13fex");
  auto h13gspAS = (TH1F *)fAS->Get("h13gsp");
  auto f2AS = new TFile("flavorProcesshists_dataAS.root");
  auto h12alldataAS = (TH1F *)f2AS->Get("h12all");
  auto h13alldataAS = (TH1F *)f2AS->Get("h13all");

  cout<<"chi2 "<<h12alldataAS->Chi2Test(h12allmcAS,"UU,NORM")<<endl;
  cout<<"KS   "<<h12alldataAS->KolmogorovTest(h12allmcAS)<<endl;

  N13allmcNS = h13allmcNS->Integral();
  N13fcrNS = h13fcrNS->Integral();
  N13fexNS = h13fexNS->Integral();
  N13gspNS = h13gspNS->Integral();
  //N12alldataNS = h12alldataNS->Integral();
  N13alldataNS = h13alldataNS->Integral();

  N13allmcAS = h13allmcAS->Integral();
  N12fcrAS = h12fcrAS->Integral();
  N12fexAS = h12fexAS->Integral();
  N12gspAS = h12gspAS->Integral();
  N13fcrAS = h13fcrAS->Integral();
  N13fexAS = h13fexAS->Integral();
  N13gspAS = h13gspAS->Integral();
  N12alldataAS = h12alldataAS->Integral();
  N13alldataAS = h13alldataAS->Integral();

  N12allmcAS = h12allmcAS->Integral();

  cout<<"N13allmcNS = "<<N13allmcNS<<endl;
  cout<<"N13fcrNS = "<<N13fcrNS<<endl;
  cout<<"N13fexNS = "<<N13fexNS<<endl;
  cout<<"N13gspNS = "<<N13gspNS<<endl;
  // cout<<"N12alldataNS = "<<N12alldataNS<<endl;
  cout<<"N13alldataNS = "<<N13alldataNS<<endl;
  cout<<"N13allmcAS = "<<N13allmcAS<<endl;
  cout<<"N13fcrAS = "<<N13fcrAS<<endl;
  cout<<"N13fexAS = "<<N13fexAS<<endl;
  cout<<"N13gspAS = "<<N13gspAS<<endl;
  cout<<"N12alldataAS = "<<N12alldataAS<<endl;
  cout<<"N13alldataAS = "<<N13alldataAS<<endl;

}



float e(float beta, float gamma)
{
  float obs1 = N12alldataAS;
  float obs2 = N13alldataAS;

  float exp1 = N13alldataAS*(N12fcrAS+beta*N12fexAS+gamma*N12gspAS)/(N13fcrAS+beta*N13fexAS+gamma*N13gspAS);
  float exp2 = N13alldataNS*(N13fcrAS+beta*N13fexAS+gamma*N13gspAS)/(N13fcrNS+beta*N13fexNS+gamma*N13gspNS);

  return (obs1-exp1)*(obs1-exp1)/abs(exp1) + (obs2-exp2)*(obs2-exp2)/abs(exp2);

  // float x = N12alldataAS/N13alldataAS - (N12fcrAS+beta*N12fexAS+gamma*N12gspAS)/(N13fcrAS+beta*N13fexAS+gamma*N13gspAS);
  // float y = N13alldataAS/N13alldataNS - (N13fcrAS+beta*N13fexAS+gamma*N13gspAS)/(N13fcrNS+beta*N13fexNS+gamma*N13gspNS);
  // return x*x+y*y;
}

void check(bool verbose, float &minimumb, float &minimumg)
{
  int N=4000;

  float betamax = 2;//0.4;//2;
  float gammamax = 2;//2;

  float betamin = -2;//-0.2;
  float gammamin = -2;//0.4;

  float minb = betamin, ming = gammamin, chi2min = e(minb, ming);
  float min = e(minb,ming);

  float bstep = (betamax-betamin)/N;
  float gstep = (gammamax-gammamin)/N;

  TH2F * h2;
  if (verbose) h2 = new TH2F("minmap","minmap",N,betamin,betamax,N,gammamin,gammamax);
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {
      float beta=betamin+bstep*i;
      float gamma=betamin+gstep*j;
      float er = e(beta,gamma);
      if (er<min) {
        min = er;
        minb = beta;
        ming = gamma;
        chi2min = er;
      }
     if (verbose)  h2->SetBinContent(i,j,er);
    }

    cout<<"minimum : "<<min<<" at beta = "<<minb<<", gamma = "<<ming<<endl;

    if (verbose) {

      float chi2e = chi2min+1;
      cout<<"chi2e"<<chi2e<<endl;
      float tol = 0.1;
      for (int i=0;i<N;i++)
      	for (int j=0;j<N;j++) {
      		if (abs(h2->GetBinContent(i,j)-chi2e)<tol)
      			h2->SetBinContent(i,j,1000);

      	}



      auto c = getc();
      h2->GetXaxis()->SetTitle("#beta");
      h2->GetYaxis()->SetTitle("#gamma");
      h2->GetYaxis()->SetRangeUser(0.9,1.5);
      h2->GetXaxis()->SetRangeUser(-0.05,0.15);
      h2->Draw("colz");


      c->SetLogz();
      c->SaveAs("plots/minmap.pdf");
    }
    minimumb = minb;
    minimumg = ming;

}

void checkerr()
{
  int n1 = N12alldataAS;
  int n2 = N13alldataAS;
  int n3 = N13alldataNS;

  TRandom *rand1 = new TRandom();
  TRandom *rand2 = new TRandom();
  TRandom *rand3 = new TRandom();

  float mb = 0, mg = 0;


  TH1F *hmb = new TH1F("hmb","hmb",100,-2,2);
  TH1F *hmg = new TH1F("hmg","hmg",100,-2,2);

  for (int i=0;i<1000;i++) {
    N12alldataAS = rand1->Poisson(n1);
    N13alldataAS = rand2->Poisson(n2);
    N13alldataNS = rand3->Poisson(n3);
    check(false, mb,mg);
    hmb->Fill(mb);
    hmg->Fill(mg);
  }

  cout<<"Mean beta  = "<<hmb->GetMean()<<"±"<<hmb->GetMeanError()<<endl;
  cout<<"Mean gamma = "<<hmg->GetMean()<<"±"<<hmg->GetMeanError()<<endl;
  Draw({hmb, hmg});

}

void PrintTable()
{

  float s1 = (N12alldataAS+N13alldataAS+N13alldataNS)/100;
  float s2 = (N12allmcAS+N13allmcAS+N13allmcNS)/100;

  float s12as = (N12fcrAS+N12gspAS+N12fexAS)/100;
  float s13as = (N13fcrAS+N13gspAS+N13fexAS)/100;
  float s13ns = (N13fcrNS+N13gspNS+N13fexNS)/100;

  cout<<setprecision(2);

  cout<<"Pair & Data & MC \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"1-2 AS & "<<N12alldataAS/s1<<"\\% & "<<N12allmcAS/s2<<"\\% \\\\"<<endl;
  cout<<"1-3 AS & "<<N13alldataAS/s1<<"\\% & "<<N13allmcAS/s2<<"\\% \\\\"<<endl;
  cout<<"1-3 NS & "<<N13alldataNS/s1<<"\\% & "<<N13allmcNS/s2<<"\\% \\\\"<<endl;
  cout<<endl;

  cout<<"Pair & FCR & GSP & FEX \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"1-2 AS & "<<N12fcrAS/s12as<<"\\% & "<<N12gspAS/s12as<<"\\% & "<<N12fexAS/s12as<<"\\% \\\\"<<endl;
  cout<<"1-3 AS & "<<N13fcrAS/s13as<<"\\% & "<<N13gspAS/s13as<<"\\% & "<<N13fexAS/s13as<<"\\% \\\\"<<endl;
  cout<<"1-3 NS & "<<N13fcrNS/s13ns<<"\\% & "<<N13gspNS/s13ns<<"\\% & "<<N13fexNS/s13ns<<"\\% \\\\"<<endl;
  cout<<endl;

}


void checknorms()
{
  //loadnumbers();
	//checknormsbeta();
	//checknormsgamma();

  loadnumbers2();
  float mb=0,mg=0;
  check(true,mb,mg);
  //checkerr();

  PrintTable();
}
