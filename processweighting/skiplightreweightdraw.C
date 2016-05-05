#include "../helpers/plotting.h"


void skiplightreweightdraw()
{
	macro m("skiplightreweightdraw_12030");

//0.266, gamma = 1.232

//0.04;//
// 1.18;//

//0.266, 0.932 - best fit for SL

	float beta  = 0.04;//0.266;//120/30 = 0.045;//110/50 = 0.28;//FEXGSP = 0.4
	float gamma = 1.2;//1.232;//120/30 = 1.182;//110/50 = 1.24;//FEXGSP = 0.4

	auto f= new TFile("skiplightreweight.root");
	auto h12all = (TH1F *)f->Get("h12all");
	auto h12fcr = (TH1F *)f->Get("h12fcr");
	auto h12fex = (TH1F *)f->Get("h12fex");
	auto h12gsp = (TH1F *)f->Get("h12gsp");
	auto hSLall = (TH1F *)f->Get("hSLall");
	auto hSLfcr = (TH1F *)f->Get("hSLfcr");
	auto hSLfex = (TH1F *)f->Get("hSLfex");
	auto hSLgsp = (TH1F *)f->Get("hSLgsp");
	auto h12data = (TH1F *)f->Get("h12data");
	auto hSLdata = (TH1F *)f->Get("hSLdata");

	auto h12dphiall = (TH1F *)f->Get("h12dphiall");
	auto h12dphifcr = (TH1F *)f->Get("h12dphifcr");
	auto h12dphifex = (TH1F *)f->Get("h12dphifex");
	auto h12dphigsp = (TH1F *)f->Get("h12dphigsp");
	auto hSLdphiall = (TH1F *)f->Get("hSLdphiall");
	auto hSLdphifcr = (TH1F *)f->Get("hSLdphifcr");
	auto hSLdphifex = (TH1F *)f->Get("hSLdphifex");
	auto hSLdphigsp = (TH1F *)f->Get("hSLdphigsp");
	auto h12dphidata = (TH1F *)f->Get("h12dphidata");
	auto hSLdphidata = (TH1F *)f->Get("hSLdphidata");

	auto h12dphiNSall = (TH1F *)f->Get("h12dphiNSall");
	auto h12dphiNSfcr = (TH1F *)f->Get("h12dphiNSfcr");
	auto h12dphiNSfex = (TH1F *)f->Get("h12dphiNSfex");
	auto h12dphiNSgsp = (TH1F *)f->Get("h12dphiNSgsp");
	auto hSLdphiNSall = (TH1F *)f->Get("hSLdphiNSall");
	auto hSLdphiNSfcr = (TH1F *)f->Get("hSLdphiNSfcr");
	auto hSLdphiNSfex = (TH1F *)f->Get("hSLdphiNSfex");
	auto hSLdphiNSgsp = (TH1F *)f->Get("hSLdphiNSgsp");
	auto h12dphiNSdata = (TH1F *)f->Get("h12dphiNSdata");
	auto hSLdphiNSdata = (TH1F *)f->Get("hSLdphiNSdata");


	auto h12ordall = (TH1F *)f->Get("h12ordall");
	auto h12ordfcr = (TH1F *)f->Get("h12ordfcr");
	auto h12ordfex = (TH1F *)f->Get("h12ordfex");
	auto h12ordgsp = (TH1F *)f->Get("h12ordgsp");
	auto hSLordall = (TH1F *)f->Get("hSLordall");
	auto hSLordfcr = (TH1F *)f->Get("hSLordfcr");
	auto hSLordfex = (TH1F *)f->Get("hSLordfex");
	auto hSLordgsp = (TH1F *)f->Get("hSLordgsp");
	auto h12orddata = (TH1F *)f->Get("h12orddata");
	auto hSLorddata = (TH1F *)f->Get("hSLorddata");


	auto h12reweighted = (TH1F *)h12all->Clone("h12reweighted");
	auto hSLreweighted = (TH1F *)hSLall->Clone("hSLreweighted");
	h12reweighted->Reset();h12reweighted->SetTitle("h12reweighted");
	hSLreweighted->Reset();hSLreweighted->SetTitle("hSLreweighted");

	auto h12dphiNSreweighted = (TH1F *)h12dphiNSall->Clone("h12dphiNSreweighted");
	auto hSLdphiNSreweighted = (TH1F *)hSLdphiNSall->Clone("hSLdphiNSreweighted");
	h12dphiNSreweighted->Reset();h12dphiNSreweighted->SetTitle("h12dphiNSreweighted");
	hSLdphiNSreweighted->Reset();hSLdphiNSreweighted->SetTitle("hSLdphiNSreweighted");

	auto h12dphireweighted = (TH1F *)h12dphiall->Clone("h12dphireweighted");
	auto hSLdphireweighted = (TH1F *)hSLdphiall->Clone("hSLdphireweighted");
	h12dphireweighted->Reset();h12dphireweighted->SetTitle("h12dphireweighted");
	hSLdphireweighted->Reset();hSLdphireweighted->SetTitle("hSLdphireweighted");

	auto h12ordreweighted = (TH1F *)h12ordall->Clone("h12ordreweighted");
	auto hSLordreweighted = (TH1F *)hSLordall->Clone("hSLordreweighted");
	h12ordreweighted->Reset();h12ordreweighted->SetTitle("h12ordreweighted");
	hSLordreweighted->Reset();hSLordreweighted->SetTitle("hSLordreweighted");


	h12reweighted->Add(h12fex,h12gsp,beta,gamma);
	h12reweighted->Add(h12fcr);

	hSLreweighted->Add(hSLfex,hSLgsp,beta,gamma);
	hSLreweighted->Add(hSLfcr);

	h12dphireweighted->Add(h12dphifex,h12dphigsp,beta,gamma);
	h12dphireweighted->Add(h12dphifcr);

	hSLdphireweighted->Add(hSLdphifex,hSLdphigsp,beta,gamma);
	hSLdphireweighted->Add(hSLdphifcr);

	h12dphiNSreweighted->Add(h12dphiNSfex,h12dphiNSgsp,beta,gamma);
	h12dphiNSreweighted->Add(h12dphiNSfcr);

	hSLdphiNSreweighted->Add(hSLdphiNSfex,hSLdphiNSgsp,beta,gamma);
	hSLdphiNSreweighted->Add(hSLdphiNSfcr);

	h12ordreweighted->Add(h12ordfex,h12ordgsp,beta,gamma);
	h12ordreweighted->Add(h12ordfcr);

	hSLordreweighted->Add(hSLordfex,hSLordgsp,beta,gamma);
	hSLordreweighted->Add(hSLordfcr);



	Normalize({h12data,h12all,hSLdata,hSLall,
			h12dphiall,hSLdphiall,h12dphidata,hSLdphidata,
			h12dphiNSall,hSLdphiNSall,h12dphiNSdata,hSLdphiNSdata,
			h12ordall,hSLordall,h12orddata,hSLorddata});
	SetMC({h12all,hSLall,h12dphiall,hSLdphiall});

	ploteffectiveentries = false;
	plotymax = 0.12;
	plotylog = false;
	plotputmean = true;
	aktstring+="PF Jets R=0.4";
  	plotsecondline = "p_{T,1}>100GeV,p_{T,2}>40GeV";
  	plotthirdline  = "#Delta#phi>2#pi/3, CSV>0.9";
  	
  	
  	plotfilenameend = "12";
  	DrawCompare(h12data,h12all);
	plotfilenameend = "SL";
  	DrawCompare(hSLdata,hSLall);

	Normalize({h12reweighted,hSLreweighted});
	SetMC({h12reweighted,hSLreweighted});

	plotfilenameend = "12";
	DrawCompare(h12data,h12reweighted);
	plotfilenameend = "SL";
	DrawCompare(hSLdata,hSLreweighted);

	plotfilenameend = "";

	// vector<int> colors = {TColor::GetColor("#FA7200"),
	// 					  TColor::GetColor("#9ACD32"),
	// 					  TColor::GetColor("#0D98BA")};//{kGreen-9, kOrange+1,kBlue+3};

	vector<int> colors = {TColor::GetColor(25,87,5),//18,58,5),
						  TColor::GetColor(255,109,24),//234,117,1),//255,117,24),
						  TColor::GetColor(77,135,232)};//{kGreen-9, kOrange+1,kBlue+3};


	auto h12stack = stackhists({h12fcr,h12fex,h12gsp},colors,"h12stack","(P)");
	DrawCompare(h12data,h12stack,"x_{J}");


	h12fex->Scale(beta); h12gsp->Scale(gamma);
	auto h12stack2 = stackhists({h12fcr,h12fex,h12gsp},colors,"h12stack2","(W)");
	DrawCompare(h12data,h12stack2,"x_{J}");


	auto hSLstack = stackhists({hSLfcr,hSLfex,hSLgsp},colors,"hSLstack","(P)");
	DrawCompare(hSLdata,hSLstack,"x_{J}");


	hSLfex->Scale(beta); hSLgsp->Scale(gamma);
	auto hSLstack2 = stackhists({hSLfcr,hSLfex,hSLgsp},colors,"hSLstack2","(W)");
	DrawCompare(hSLdata,hSLstack2,"x_{J}");





	plotputmean = false;
	plotylog = false;
	plotymax = 1.1;
	plotymin = 1E-5;


	Normalize({hSLordreweighted,h12ordreweighted});
	SetMC({hSLordreweighted,h12ordreweighted});

	DrawCompare(h12orddata,h12ordreweighted);
	DrawCompare(hSLorddata,hSLordreweighted);

	auto h12ordstack = stackhists({h12ordfcr,h12ordfex,h12ordgsp},colors,"h12ordstack","(P)");
	DrawCompare(h12orddata,h12ordstack,"skiplight order");


	h12ordfex->Scale(beta); h12ordgsp->Scale(gamma);
	auto h12ordstack2 = stackhists({h12ordfcr,h12ordfex,h12ordgsp},colors,"h12ordstack2","(W)");
	DrawCompare(h12orddata,h12ordstack2,"skiplight order");


	auto hSLordstack = stackhists({hSLordfcr,hSLordfex,hSLordgsp},colors,"hSLordstack","(P)");
	DrawCompare(hSLorddata,hSLordstack,"skiplight order");


	hSLordfex->Scale(beta); hSLordgsp->Scale(gamma);
	auto hSLordstack2 = stackhists({hSLordfcr,hSLordfex,hSLordgsp},colors,"hSLordstack2","(W)");
	DrawCompare(hSLorddata,hSLordstack2,"skiplight order");






	plotylog = true;
	plotymin = 9999;
	plotputmean = false;
	plotputwidth = true;
	plotymax = 0.5;
	plotthirdline  = "CSV>0.9";

	DrawCompare(h12dphidata,h12dphiall,"#Delta#phi");
	DrawCompare(hSLdphidata,hSLdphiall,"#Delta#phi");


	Normalize({hSLdphireweighted,h12dphireweighted});
	SetMC({hSLdphireweighted,h12dphireweighted});

	DrawCompare(h12dphidata,h12dphireweighted,"#Delta#phi");
	DrawCompare(hSLdphidata,hSLdphireweighted,"#Delta#phi");

	//plotymax = 0.5;

	auto h12dphistack = stackhists({h12dphifcr,h12dphifex,h12dphigsp},colors,"h12dphistack","(P)");
	DrawCompare(h12dphidata,h12dphistack,"#Delta#phi");

	h12dphifex->Scale(beta); h12dphigsp->Scale(gamma);
	auto h12dphistack2 = stackhists({h12dphifcr,h12dphifex,h12dphigsp},colors,"h12dphistack2","(W)");
	DrawCompare(h12dphidata,h12dphistack2,"#Delta#phi");


	auto hSLdphistack = stackhists({hSLdphifcr,hSLdphifex,hSLdphigsp},colors,"hSLdphistack","(P)");
	DrawCompare(hSLdphidata,hSLdphistack,"#Delta#phi");

	hSLdphifex->Scale(beta); hSLdphigsp->Scale(gamma);
	auto hSLdphistack2 = stackhists({hSLdphifcr,hSLdphifex,hSLdphigsp},colors,"hSLdphistack2","(W)");
	DrawCompare(hSLdphidata,hSLdphistack2,"#Delta#phi");



	plotylog = false;
	plotymax = 0.25;
//only interesting dphi region
	DrawCompare(h12dphiNSdata,h12dphiNSall);
	DrawCompare(hSLdphiNSdata,hSLdphiNSall);


	Normalize({hSLdphiNSreweighted,h12dphiNSreweighted});
	SetMC({hSLdphiNSreweighted,h12dphiNSreweighted});

	DrawCompare(h12dphiNSdata,h12dphiNSreweighted);
	DrawCompare(hSLdphiNSdata,hSLdphiNSreweighted);



	auto h12dphiNSstack = stackhists({h12dphiNSfcr,h12dphiNSfex,h12dphiNSgsp},colors,"h12dphiNSstack","(P)");
	DrawCompare(h12dphiNSdata,h12dphiNSstack,"#Delta#phi");

	h12dphiNSfex->Scale(beta); h12dphiNSgsp->Scale(gamma);
	auto h12dphiNSstack2 = stackhists({h12dphiNSfcr,h12dphiNSfex,h12dphiNSgsp},colors,"h12dphiNSstack2","(W)");
	DrawCompare(h12dphiNSdata,h12dphiNSstack2,"#Delta#phi");


	auto hSLdphiNSstack = stackhists({hSLdphiNSfcr,hSLdphiNSfex,hSLdphiNSgsp},colors,"hSLdphiNSstack","(P)");
	DrawCompare(hSLdphiNSdata,hSLdphiNSstack,"#Delta#phi");

	hSLdphiNSfex->Scale(beta); hSLdphiNSgsp->Scale(gamma);
	auto hSLdphiNSstack2 = stackhists({hSLdphiNSfcr,hSLdphiNSfex,hSLdphiNSgsp},colors,"hSLdphiNSstack2","(W)");
	DrawCompare(hSLdphiNSdata,hSLdphiNSstack2,"#Delta#phi");



// 	int nearsidebin1 = 1;
// 	int nearsidebin2 = hSLdphifex->GetNbinsX()/3;

// 	int awaysidebin1 = hSLdphifex->GetNbinsX()/3*2;
// 	int awaysidebin2 = hSLdphifex->GetNbinsX()-1;

// 	float NnsFEX  = hSLdphifex->Integral(nearsidebin1,nearsidebin2);
// 	float NnsGSP  = hSLdphigsp->Integral(nearsidebin1,nearsidebin2);
// 	float NasFEX  = hSLdphifex->Integral(awaysidebin1,awaysidebin2);
// 	float NasGSP  = hSLdphigsp->Integral(awaysidebin1,awaysidebin2);

// 	float NnsData = hSLdphidata->Integral(nearsidebin1,nearsidebin2);
// 	float NasData = hSLdphidata->Integral(awaysidebin1,awaysidebin2);

// 	float alpha = (NnsData - NnsGSP*(1+NasFEX/NasGSP))/(NnsFEX - NnsGSP*NasFEX/NasGSP); 
// 	float beta2 = ((1-alpha)*NasFEX+NasGSP)/NasGSP;

// if (alpha<0) alpha =0;

// 	cout<<"!!!"<<alpha<<" "<<beta2<<endl;



// 	hSLdphifex->Scale(alpha); hSLdphigsp->Scale(beta2);
// 	auto hSLdphistack2_new = stackhists({hSLdphifcr,hSLdphifex,hSLdphigsp},colors,"hSLdphistack2_new");
// 	DrawCompare(hSLdphidata,hSLdphistack2_new,"#Delta#phi");


// 	cout<<" NasFEX+NasGSP = "<<NasFEX+NasGSP<<" ?= "<<
// 		hSLdphifex->Integral(awaysidebin1,awaysidebin2)+hSLdphigsp->Integral(awaysidebin1,awaysidebin2)<<endl;


// 	ploteffectiveentries = false;
// 	plotymax = 0.12;
// 	plotylog = false;
// 	plotputmean = true;
// 	aktstring+="PF Jets R=0.4";
//   	plotsecondline = "p_{T,1}>120GeV,p_{T,2}>30GeV";
//   	plotthirdline  = "#Delta#phi>2#pi/3, CSV>0.9";






// 	auto h12reweightedtwice = (TH1F *)h12all->Clone("h12reweightedtwice");
// 	auto hSLreweightedtwice = (TH1F *)hSLall->Clone("hSLreweightedtwice");
// 	h12reweightedtwice->Reset();h12reweightedtwice->SetTitle("h12reweightedtwice");
// 	hSLreweightedtwice->Reset();hSLreweightedtwice->SetTitle("hSLreweightedtwice");

// 	h12reweightedtwice->Add(h12fex,h12gsp,alpha,beta2);
// 	h12reweightedtwice->Add(h12fcr);
// 	cout<<"h12reweightedtwice "<<h12reweightedtwice->GetMean()<<endl;

// 	hSLreweightedtwice->Add(hSLfex,hSLgsp,alpha,beta2);
// 	hSLreweightedtwice->Add(hSLfcr);
// 	cout<<"hSLreweightedtwice "<<hSLreweightedtwice->GetMean()<<endl;

// 	Normalize({h12reweightedtwice,hSLreweightedtwice});
// 	SetMC({h12reweightedtwice,hSLreweightedtwice});

// 	plotfilenameend = "12";
// 	DrawCompare(h12data,h12reweightedtwice);
// 	plotfilenameend = "SL";
// 	DrawCompare(hSLdata,hSLreweightedtwice);
// 	plotfilenameend="";



// 	hSLfex->Scale(alpha); hSLgsp->Scale(beta2);
// 	auto hSLstack2_new = stackhists({hSLfcr,hSLfex,hSLgsp},colors,"hSLstack2_new");
// 	DrawCompare(hSLdata,hSLstack2_new,"x_{J}");

// 	h12fex->Scale(alpha); h12gsp->Scale(beta2);
// 	auto h12stack2_new = stackhists({h12fcr,h12fex,h12gsp},colors,"h12stack2_new");
// 	DrawCompare(h12data,h12stack2_new,"x_{J}","___");






}