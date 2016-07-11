#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/physics.h"
#include "../helpers/config.h"

bool dphiAScut = true;

bool dphicut(float dphi)
{
  float pi = acos(-1);
  float pi13 = 1./3*pi;
  float pi23 = 2./3*pi;
  
  return dphiAScut ? dphi>pi23 : dphi<pi13;

}

float jtpt1min = 100;//120;
float jtpt2min = 40;//30;

void flavorProcessFile(bool mc, TString filename, TString filenameout)
{
	auto f = new TFile(filename);
	//auto nt=(TTree *)f->Get("nt");

	auto f2 = new TFile(filenameout ,"recreate");

	auto h12all = geth("h12all","ALL jet12");
	auto h12fcr = geth("h12fcr","FCR jet12");
	auto h12fex = geth("h12fex","FEX jet12");
	auto h12gsp = geth("h12gsp","GSP jet12");
	
	auto h13all = geth("h13all","ALL jet13");
	auto h13fcr = geth("h13fcr","FCR jet13");
	auto h13fex = geth("h13fex","FEX jet13");
	auto h13gsp = geth("h13gsp","GSP jet13");
		
	// auto h23all = geth("h23all","ALL jet23");
	// auto h23fcr = geth("h23fcr","FCR jet23");
	// auto h23fex = geth("h23fex","FEX jet23");
	// auto h23gsp = geth("h23gsp","GSP jet23");


	unordered_set<int> eventstomiss = {646529,664401,720507};

	Fill(f, [&] (dict &m) {
		float w = m["weight"];
		if (!mc && w<0) return;
		if (mc && m["pthat"]<pthatcut) return;
		if (mc && eventstomiss.find((int)m["event"])!=eventstomiss.end()) return;

		if (m["jtpt1"]>jtpt1min && m["jtpt2"]>jtpt2min && m["jtpt3"]>jtpt2min && dphicut(m["dphi21"]) && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_2"]>0.9 && !(m["discr_csvV1_3"]>0.5)) {
			float xj = m["jtpt2"]/m["jtpt1"];
			h12all->Fill(xj,w);

			if (mc && m["bProdCode"]==1) h12fcr->Fill(xj,w);
			if (mc && m["bProdCode"]==2) h12fex->Fill(xj,w);
			if (mc && m["bProdCode"]==0) h12gsp->Fill(xj,w);
		}

		if (m["jtpt1"]>jtpt1min && m["jtpt3"]>jtpt2min && m["jtpt2"]>jtpt2min && dphicut(m["dphi31"]) && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_3"]>0.9 && !(m["discr_csvV1_2"]>0.5)) {
			float xj = m["jtpt3"]/m["jtpt1"];
			h13all->Fill(xj,w);

			if (mc && m["bProdCode"]==1) h13fcr->Fill(xj,w);
			if (mc && m["bProdCode"]==2) h13fex->Fill(xj,w);
			if (mc && m["bProdCode"]==0) h13gsp->Fill(xj,w);
		}

		// if (m["jtpt1"]>jtpt1min && m["jtpt3"]>jtpt2min && dphicut(m["dphi31"]) && m["discr_csvV1_1"]<0.5 && m["discr_csvV1_2"]>0.9 && m["discr_csvV1_3"]>0.9) {
		// 	float xj = m["jtpt3"]/m["jtpt1"];
		// 	h23all->Fill(xj,w);

		// 	if (mc && m["bProdCode"]==1) h23fcr->Fill(xj,w);
		// 	if (mc && m["bProdCode"]==2) h23fex->Fill(xj,w);
		// 	if (mc && m["bProdCode"]==0) h23gsp->Fill(xj,w);
		// }

	});

	f2->cd();

	h12all->Write();
	h12fcr->Write();
	h12fex->Write();
	h12gsp->Write();
	
	h13all->Write();
	h13fcr->Write();
	h13fex->Write();
	h13gsp->Write();
	
	// h23all->Write();
	// h23fcr->Write();
	// h23fex->Write();
	// h23gsp->Write();

	f2->Close();
	f->Close();



}

void flavorProcess()
{
	dphiAScut = true;
	flavorProcessFile(true,config.getFileName_djt("mcppbfa"), "flavorProcesshists_mcAS.root");
	flavorProcessFile(false,config.getFileName_djt("dtppjpf"), "flavorProcesshists_dataAS.root");

	dphiAScut = false;
	flavorProcessFile(true,config.getFileName_djt("mcppbfa"), "flavorProcesshists_mcNS.root");
	flavorProcessFile(false,config.getFileName_djt("dtppjpf"), "flavorProcesshists_dataNS.root");

}
