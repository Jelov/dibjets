#include "../helpers/plotting.h"
#include "../helpers/looptuple.h"
#include "../helpers/config.h"
#include "../helpers/physics.h"

bool pthatcondition(float pthat)
{
  return pthat>pthatcut;//pthat==120;
	//return true;
}

float jtpt1min = 100;//120;
float jtpt2min = 40;//30;

void skiplightreweight()
{
	TString outfilename = "skiplightreweight.root";

	auto fmc = config.getfile_djt("mcppbfa");
	auto fdt = config.getfile_djt("dtppjpf");

	auto fout = new TFile(outfilename ,"recreate");


	auto h12all = geth("h12all","ALL"); // jet12
	auto h12fcr = geth("h12fcr","FCR"); // jet12
	auto h12fex = geth("h12fex","FEX"); // jet12
	auto h12gsp = geth("h12gsp","GSP"); // jet12
	
	auto hSLall = geth("hSLall","ALL"); // jetSL
	auto hSLfcr = geth("hSLfcr","FCR"); // jetSL
	auto hSLfex = geth("hSLfex","FEX"); // jetSL
	auto hSLgsp = geth("hSLgsp","GSP"); // jetSL

	auto h12data = geth("h12data","DATA"); // jet12
	auto hSLdata = geth("hSLdata","DATA"); // jetSL

	buildxmin = 0;
	buildxmax = 3.142;

	auto h12dphiall = geth("h12dphiall","ALL"); // dphi jet12
	auto h12dphifcr = geth("h12dphifcr","FCR"); // dphi jet12
	auto h12dphifex = geth("h12dphifex","FEX"); // dphi jet12
	auto h12dphigsp = geth("h12dphigsp","GSP"); // dphi jet12
	auto hSLdphiall = geth("hSLdphiall","ALL"); // dphi jetSL
	auto hSLdphifcr = geth("hSLdphifcr","FCR"); // dphi jetSL
	auto hSLdphifex = geth("hSLdphifex","FEX"); // dphi jetSL
	auto hSLdphigsp = geth("hSLdphigsp","GSP"); // dphi jetSL

	auto h12dphidata = geth("h12dphidata","DATA"); // jetSL dphi
	auto hSLdphidata = geth("hSLdphidata","DATA"); // jetSL dphi

	buildxmin = 2./3*3.1416;
	buildxmax = 3.1416;

	auto h12dphiNSall = geth("h12dphiNSall","ALL"); // dphi jet12
	auto h12dphiNSfcr = geth("h12dphiNSfcr","FCR"); // dphi jet12
	auto h12dphiNSfex = geth("h12dphiNSfex","FEX"); // dphi jet12
	auto h12dphiNSgsp = geth("h12dphiNSgsp","GSP"); // dphi jet12
	auto hSLdphiNSall = geth("hSLdphiNSall","ALL"); // dphi jetSL
	auto hSLdphiNSfcr = geth("hSLdphiNSfcr","FCR"); // dphi jetSL
	auto hSLdphiNSfex = geth("hSLdphiNSfex","FEX"); // dphi jetSL
	auto hSLdphiNSgsp = geth("hSLdphiNSgsp","GSP"); // dphi jetSL

	auto h12dphiNSdata = geth("h12dphiNSdata","DATA"); // jetSL dphi
	auto hSLdphiNSdata = geth("hSLdphiNSdata","DATA"); // jetSL dphi

	buildnbins = 6;
	buildxmin = 1;
	buildxmax = 7;

	auto h12ordall = geth("h12ordall","ALL"); // ord jet12
	auto h12ordfcr = geth("h12ordfcr","FCR"); // ord jet12
	auto h12ordfex = geth("h12ordfex","FEX"); // ord jet12
	auto h12ordgsp = geth("h12ordgsp","GSP"); // ord jet12
	auto hSLordall = geth("hSLordall","ALL"); // ord jetSL
	auto hSLordfcr = geth("hSLordfcr","FCR"); // ord jetSL
	auto hSLordfex = geth("hSLordfex","FEX"); // ord jetSL
	auto hSLordgsp = geth("hSLordgsp","GSP"); // ord jetSL
	auto h12orddata = geth("h12orddata","DATA"); // jetSL ord
	auto hSLorddata = geth("hSLorddata","DATA"); // jetSL ord

	
	float pi23 = 2./3*acos(-1);

	vector<TString> mcvariables = {"jtpt1","jtpt2","jtptSL","discr_csvV1_1","discr_csvV1_2","discr_csvV1_SL","SLord",
	"dphi21","dphiSL1","weight","bProdCode","pthat"};
	vector<TString> dtvariables = {"jtpt1","jtpt2","jtptSL","discr_csvV1_1","discr_csvV1_2","discr_csvV1_SL","SLord",
	"dphi21","dphiSL1","weight","hltPFJet80"};

	Fill(fmc, mcvariables, [&] (dict &m) -> void {
		float w = m["weight"];

		if (m["jtpt1"]>jtpt1min && m["jtpt2"]>jtpt2min && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_2"]>0.9) {
			float dphi21= m["dphi21"];

			if ( pthatcondition(m["pthat"])) { h12dphiall->Fill(dphi21,w); h12dphiNSall->Fill(dphi21,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==1) { h12dphifcr->Fill(dphi21,w); h12dphiNSfcr->Fill(dphi21,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==2) { h12dphifex->Fill(dphi21,w); h12dphiNSfex->Fill(dphi21,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==0) { h12dphigsp->Fill(dphi21,w); h12dphiNSgsp->Fill(dphi21,w);}

			if (dphi21>pi23) {
				float xj = m["jtpt2"]/m["jtpt1"];
				float slord = m["SLord"];
				if ( pthatcondition(m["pthat"])) {h12all->Fill(xj,w);h12ordall->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==1) {h12fcr->Fill(xj,w);h12ordfcr->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==2) {h12fex->Fill(xj,w);h12ordfex->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==0) {h12gsp->Fill(xj,w);h12ordgsp->Fill(slord,w);}
			}
		}
//well... SL is tagged by definition...
		if (m["jtpt1"]>jtpt1min && m["jtptSL"]>jtpt2min && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_SL"]>0.9) {
			float dphiSL = m["dphiSL1"];

			if ( pthatcondition(m["pthat"])) { hSLdphiall->Fill(dphiSL,w); hSLdphiNSall->Fill(dphiSL,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==1) { hSLdphifcr->Fill(dphiSL,w); hSLdphiNSfcr->Fill(dphiSL,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==2) { hSLdphifex->Fill(dphiSL,w); hSLdphiNSfex->Fill(dphiSL,w);}
			if ( pthatcondition(m["pthat"]) && m["bProdCode"]==0) { hSLdphigsp->Fill(dphiSL,w); hSLdphiNSgsp->Fill(dphiSL,w);}

			if (dphiSL>pi23) {
				float xj = m["jtptSL"]/m["jtpt1"];
				float slord = m["SLord"];
				if ( pthatcondition(m["pthat"])) {hSLall->Fill(xj,w); hSLordall->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==1) {hSLfcr->Fill(xj,w);hSLordfcr->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==2) {hSLfex->Fill(xj,w);hSLordfex->Fill(slord,w);}
				if ( pthatcondition(m["pthat"]) && m["bProdCode"]==0) {hSLgsp->Fill(xj,w);hSLordgsp->Fill(slord,w);}
			}
		}

	});

	Fill(fdt, dtvariables, [&] (dict &m) -> void {
		float w = m["weight"];
		if (m["jtpt1"]>jtpt1min && m["jtpt2"]>jtpt2min && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_2"]>0.9) {
			float xj = m["jtpt2"]/m["jtpt1"];
			float slord = m["SLord"];
			float dphi21 = m["dphi21"];

			h12dphidata->Fill(dphi21,w);
			h12dphiNSdata->Fill(dphi21,w);

			if (dphi21>pi23) {
				h12data->Fill(xj,w);
				h12orddata->Fill(slord,w);
			}
		}

		if (m["jtpt1"]>jtpt1min && m["jtptSL"]>jtpt2min && m["discr_csvV1_1"]>0.9 && m["discr_csvV1_SL"]>0.9) {
			float xj = m["jtptSL"]/m["jtpt1"];
			float slord = m["SLord"];			
			float dphiSL1 = m["dphiSL1"];
			hSLdphidata->Fill(dphiSL1,w);
			hSLdphiNSdata->Fill(dphiSL1,w);

			if (dphiSL1>pi23) {
				hSLdata->Fill(xj,w);
				hSLorddata->Fill(slord,w);
			}
		}

	});

	fout->cd();
	h12all->Write();
	h12fcr->Write();
	h12fex->Write();
	h12gsp->Write();
	hSLall->Write();
	hSLfcr->Write();
	hSLfex->Write();
	hSLgsp->Write();
	h12data->Write();
	hSLdata->Write();
	h12dphiall->Write();
	h12dphifcr->Write();
	h12dphifex->Write();
	h12dphigsp->Write();
	hSLdphiall->Write();
	hSLdphifcr->Write();
	hSLdphifex->Write();
	hSLdphigsp->Write();
	h12dphidata->Write();
	hSLdphidata->Write();

	h12dphiNSall->Write(); 
	h12dphiNSfcr->Write(); 
	h12dphiNSfex->Write(); 
	h12dphiNSgsp->Write(); 
	hSLdphiNSall->Write(); 
	hSLdphiNSfcr->Write(); 
	hSLdphiNSfex->Write(); 
	hSLdphiNSgsp->Write(); 
	h12dphiNSdata->Write();
	hSLdphiNSdata->Write();

	h12ordall->Write();
	h12ordfcr->Write();
	h12ordfex->Write();
	h12ordgsp->Write();
	hSLordall->Write();
	hSLordfcr->Write();
	hSLordfex->Write();
	hSLordgsp->Write();
	h12orddata->Write();
	hSLorddata->Write();
	fdt->Close();
	fmc->Close();



}
