#include "../helpers/config.h"
float etacut = 1.5;

void iterTrigEffCorr(bool write =true)
{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Apr 18 10:15:19 2016 by ROOT version6.06/00)
//   from TTree nt/ntdj
//   found on file: ./merged_dtPbjakPu4PF_djt.root
//////////////////////////////////////////////////////////
  TH1::SetDefaultSumw2();
  
  
  TFile *f = config.getfile_djt("dtPbjcl");//new TFile("../../ntuples/eta1p5_jecv2/dtPbjclakPu4PF_djt.root");//"./merged_dtPbjakPu4PF_djt.root");
  TNtuple *nt = (TNtuple *) f->Get("nt");

//Declaration of leaves types
   Float_t         run;
   Float_t         lumi;
   Float_t         event;
   Float_t         prew;
   Float_t         triggermatched;
   Float_t         bin;
   Float_t         vz;
   Float_t         hiHF;
   Float_t         hltCSV60;
   Float_t         hltCSV80;
   Float_t         hltCaloJet40;
   Float_t         hltCaloJet60;
   Float_t         hltCaloJet80;
   Float_t         hltPFJet60;
   Float_t         hltPFJet80;
   Float_t         dijet;
   Float_t         hltCalo60jtpt;
   Float_t         hltCalo60jtphi;
   Float_t         hltCalo60jteta;
   Float_t         hltCalo80jtpt;
   Float_t         hltCalo80jtphi;
   Float_t         hltCalo80jteta;
   Float_t         hltCSV60jtpt;
   Float_t         hltCSV60jtphi;
   Float_t         hltCSV60jteta;
   Float_t         hltCSV80jtpt;
   Float_t         hltCSV80jtphi;
   Float_t         hltCSV80jteta;
   Float_t         rawpt1;
   Float_t         jtpt1;
   Float_t         jtphi1;
   Float_t         jteta1;
   Float_t         discr_csvV1_1;
   Float_t         svtxm1;
   Float_t         discr_prob1;
   Float_t         svtxdls1;
   Float_t         svtxpt1;
   Float_t         svtxntrk1;
   Float_t         nsvtx1;
   Float_t         nselIPtrk1;
   Float_t         rawpt2;
   Float_t         jtpt2;
   Float_t         jtphi2;
   Float_t         jteta2;
   Float_t         discr_csvV1_2;
   Float_t         svtxm2;
   Float_t         discr_prob2;
   Float_t         svtxdls2;
   Float_t         svtxpt2;
   Float_t         svtxntrk2;
   Float_t         nsvtx2;
   Float_t         nselIPtrk2;
   Float_t         dphi21;
   Float_t         rawpt3;
   Float_t         jtpt3;
   Float_t         jtphi3;
   Float_t         jteta3;
   Float_t         discr_csvV1_3;
   Float_t         svtxm3;
   Float_t         discr_prob3;
   Float_t         svtxdls3;
   Float_t         svtxpt3;
   Float_t         svtxntrk3;
   Float_t         nsvtx3;
   Float_t         nselIPtrk3;
   Float_t         dphi31;
   Float_t         dphi32;
   Float_t         SLord;
   Float_t         rawptSL;
   Float_t         jtptSL;
   Float_t         jtphiSL;
   Float_t         jtetaSL;
   Float_t         discr_csvV1_SL;
   Float_t         svtxmSL;
   Float_t         discr_probSL;
   Float_t         svtxdlsSL;
   Float_t         svtxptSL;
   Float_t         svtxntrkSL;
   Float_t         nsvtxSL;
   Float_t         nselIPtrkSL;
   Float_t         dphiSL1;
   Float_t         weight;

   // Set branch addresses.
//   nt->SetBranchAddress("run",&run);
//   nt->SetBranchAddress("lumi",&lumi);
//   nt->SetBranchAddress("event",&event);
//   nt->SetBranchAddress("prew",&prew);
//   nt->SetBranchAddress("triggermatched",&triggermatched);
   nt->SetBranchAddress("bin",&bin);
//   nt->SetBranchAddress("vz",&vz);
//   nt->SetBranchAddress("hiHF",&hiHF);
   nt->SetBranchAddress("hltCSV60",&hltCSV60);
   nt->SetBranchAddress("hltCSV80",&hltCSV80);
//   nt->SetBranchAddress("hltCaloJet40",&hltCaloJet40);
   nt->SetBranchAddress("hltCaloJet60",&hltCaloJet60);
   nt->SetBranchAddress("hltCaloJet80",&hltCaloJet80);
//   nt->SetBranchAddress("hltPFJet60",&hltPFJet60);
//   nt->SetBranchAddress("hltPFJet80",&hltPFJet80);
//   nt->SetBranchAddress("dijet",&dijet);
//   nt->SetBranchAddress("hltCalo60jtpt",&hltCalo60jtpt);
//   nt->SetBranchAddress("hltCalo60jtphi",&hltCalo60jtphi);
//   nt->SetBranchAddress("hltCalo60jteta",&hltCalo60jteta);
//   nt->SetBranchAddress("hltCalo80jtpt",&hltCalo80jtpt);
//   nt->SetBranchAddress("hltCalo80jtphi",&hltCalo80jtphi);
//   nt->SetBranchAddress("hltCalo80jteta",&hltCalo80jteta);
//   nt->SetBranchAddress("hltCSV60jtpt",&hltCSV60jtpt);
//   nt->SetBranchAddress("hltCSV60jtphi",&hltCSV60jtphi);
//   nt->SetBranchAddress("hltCSV60jteta",&hltCSV60jteta);
//   nt->SetBranchAddress("hltCSV80jtpt",&hltCSV80jtpt);
//   nt->SetBranchAddress("hltCSV80jtphi",&hltCSV80jtphi);
//   nt->SetBranchAddress("hltCSV80jteta",&hltCSV80jteta);
//   nt->SetBranchAddress("rawpt1",&rawpt1);
   nt->SetBranchAddress("jtpt1",&jtpt1);
//   nt->SetBranchAddress("jtphi1",&jtphi1);
   nt->SetBranchAddress("jteta1",&jteta1);
   nt->SetBranchAddress("discr_csvV1_1",&discr_csvV1_1);
//   nt->SetBranchAddress("svtxm1",&svtxm1);
//   nt->SetBranchAddress("discr_prob1",&discr_prob1);
//   nt->SetBranchAddress("svtxdls1",&svtxdls1);
//   nt->SetBranchAddress("svtxpt1",&svtxpt1);
//   nt->SetBranchAddress("svtxntrk1",&svtxntrk1);
//   nt->SetBranchAddress("nsvtx1",&nsvtx1);
//   nt->SetBranchAddress("nselIPtrk1",&nselIPtrk1);
//   nt->SetBranchAddress("rawpt2",&rawpt2);
//   nt->SetBranchAddress("jtpt2",&jtpt2);
//   nt->SetBranchAddress("jtphi2",&jtphi2);
//   nt->SetBranchAddress("jteta2",&jteta2);
//   nt->SetBranchAddress("discr_csvV1_2",&discr_csvV1_2);
//   nt->SetBranchAddress("svtxm2",&svtxm2);
//   nt->SetBranchAddress("discr_prob2",&discr_prob2);
//   nt->SetBranchAddress("svtxdls2",&svtxdls2);
//   nt->SetBranchAddress("svtxpt2",&svtxpt2);
//   nt->SetBranchAddress("svtxntrk2",&svtxntrk2);
//   nt->SetBranchAddress("nsvtx2",&nsvtx2);
//   nt->SetBranchAddress("nselIPtrk2",&nselIPtrk2);
//   nt->SetBranchAddress("dphi21",&dphi21);
//   nt->SetBranchAddress("rawpt3",&rawpt3);
//   nt->SetBranchAddress("jtpt3",&jtpt3);
//   nt->SetBranchAddress("jtphi3",&jtphi3);
//   nt->SetBranchAddress("jteta3",&jteta3);
//   nt->SetBranchAddress("discr_csvV1_3",&discr_csvV1_3);
//   nt->SetBranchAddress("svtxm3",&svtxm3);
//   nt->SetBranchAddress("discr_prob3",&discr_prob3);
//   nt->SetBranchAddress("svtxdls3",&svtxdls3);
//   nt->SetBranchAddress("svtxpt3",&svtxpt3);
//   nt->SetBranchAddress("svtxntrk3",&svtxntrk3);
//   nt->SetBranchAddress("nsvtx3",&nsvtx3);
//   nt->SetBranchAddress("nselIPtrk3",&nselIPtrk3);
//   nt->SetBranchAddress("dphi31",&dphi31);
//   nt->SetBranchAddress("dphi32",&dphi32);
//   nt->SetBranchAddress("SLord",&SLord);
//   nt->SetBranchAddress("rawptSL",&rawptSL);
//   nt->SetBranchAddress("jtptSL",&jtptSL);
//   nt->SetBranchAddress("jtphiSL",&jtphiSL);
//   nt->SetBranchAddress("jtetaSL",&jtetaSL);
//   nt->SetBranchAddress("discr_csvV1_SL",&discr_csvV1_SL);
//   nt->SetBranchAddress("svtxmSL",&svtxmSL);
//   nt->SetBranchAddress("discr_probSL",&discr_probSL);
//   nt->SetBranchAddress("svtxdlsSL",&svtxdlsSL);
//   nt->SetBranchAddress("svtxptSL",&svtxptSL);
//   nt->SetBranchAddress("svtxntrkSL",&svtxntrkSL);
//   nt->SetBranchAddress("nsvtxSL",&nsvtxSL);
//   nt->SetBranchAddress("nselIPtrkSL",&nselIPtrkSL);
//   nt->SetBranchAddress("dphiSL1",&dphiSL1);
   nt->SetBranchAddress("weight",&weight);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// nt->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = nt->GetEntries();


  TH1F *hCentJet = new TH1F("hCentJet","hCentJet",20,0,100);
  TH1F *hCentBJet = new TH1F("hCentBJet","hCentBJet",20,0,100);
  TH1F *hCentBJetEff = new TH1F("hCentBJetEff","hCentBJetEff",20,0,100);
  
  TH1F *hw1CentJet = new TH1F("hw1CentJet","hw1CentJet",20,0,100);
  TH1F *hw1CentBJet = new TH1F("hw1CentBJet","hw1CentBJet",20,0,100);
  TH1F *hw1CentBJetEff = new TH1F("hw1CentBJetEff","hw1CentBJetEff",20,0,100);

  TH1F *hw2CentJet = new TH1F("hw2CentJet","hw2CentJet",20,0,100);
  TH1F *hw2CentBJet = new TH1F("hw2CentBJet","hw2CentBJet",20,0,100);
  TH1F *hw2CentBJetEff = new TH1F("hw2CentBJetEff","hw2CentBJetEff",20,0,100);

  TH1F *hw3CentJet = new TH1F("hw3CentJet","hw3CentJet",20,0,100);
  TH1F *hw3CentBJet = new TH1F("hw3CentBJet","hw3CentBJet",20,0,100);
  TH1F *hw3CentBJetEff = new TH1F("hw3CentBJetEff","hw3CentBJetEff",20,0,100);
  
  TH1F *hEtaJet = new TH1F("hEtaJet","hEtaJet",20,-etacut,etacut);
  TH1F *hEtaBJet = new TH1F("hEtaBJet","hEtaBJet",20,-etacut,etacut);
  TH1F *hEtaBJetEff = new TH1F("hEtaBJetEff","hEtaBJetEff",20,-etacut,etacut);
  
  TH1F *hw1EtaJet = new TH1F("hw1EtaJet","hw1EtaJet",20,-etacut,etacut);
  TH1F *hw1EtaBJet = new TH1F("hw1EtaBJet","hw1EtaBJet",20,-etacut,etacut);
  TH1F *hw1EtaBJetEff = new TH1F("hw1EtaBJetEff","hw1EtaBJetEff",20,-etacut,etacut);
  
  TH1F *hw2EtaJet = new TH1F("hw2EtaJet","hw2EtaJet",20,-etacut,etacut);
  TH1F *hw2EtaBJet = new TH1F("hw2EtaBJet","hw2EtaBJet",20,-etacut,etacut);
  TH1F *hw2EtaBJetEff = new TH1F("hw2EtaBJetEff","hw2EtaBJetEff",20,-etacut,etacut);
  
  TH1F *hw3EtaJet = new TH1F("hw3EtaJet","hw3EtaJet",20,-etacut,etacut);
  TH1F *hw3EtaBJet = new TH1F("hw3EtaBJet","hw3EtaBJet",20,-etacut,etacut);
  TH1F *hw3EtaBJetEff = new TH1F("hw3EtaBJetEff","hw3EtaBJetEff",20,-etacut,etacut);
  
  
  TH1F *hPtJet = new TH1F("hPtJet","hPtJet",20,100.,200.);
  TH1F *hPtBJet = new TH1F("hPtBJet","hPtBJet",20,100.,200.);
  TH1F *hPtBJetEff = new TH1F("hPtBJetEff","hPtBJetEff",20,100.,200.);

  TH1F *hw1PtJet = new TH1F("hw1PtJet","hw1PtJet",20,100.,200.);
  TH1F *hw1PtBJet = new TH1F("hw1PtBJet","hw1PtBJet",20,100.,200.);
  TH1F *hw1PtBJetEff = new TH1F("hw1PtBJetEff","hw1PtBJetEff",20,100.,200.);

    TH1F *hw2PtJet = new TH1F("hw2PtJet","hw2PtJet",20,100.,200.);
  TH1F *hw2PtBJet = new TH1F("hw2PtBJet","hw2PtBJet",20,100.,200.);
  TH1F *hw2PtBJetEff = new TH1F("hw2PtBJetEff","hw2PtBJetEff",20,100.,200.);

    TH1F *hw3PtJet = new TH1F("hw3PtJet","hw3PtJet",20,100.,200.);
  TH1F *hw3PtBJet = new TH1F("hw3PtBJet","hw3PtBJet",20,100.,200.);
  TH1F *hw3PtBJetEff = new TH1F("hw3PtBJetEff","hw3PtBJetEff",20,100.,200.);


  TH1F *hRCentJet = new TH1F("hRCentJet","hRCentJet",10,0,100);
  TH1F *hRCentBJet = new TH1F("hRCentBJet","hRCentBJet",10,0,100);
  TH1F *hRCentBJetEff = new TH1F("hRCentBJetEff","hRCentBJetEff",10,0,100);
  
  TH1F *hRw1CentJet = new TH1F("hRw1CentJet","hRw1CentJet",10,0,100);
  TH1F *hRw1CentBJet = new TH1F("hRw1CentBJet","hRw1CentBJet",10,0,100);
  TH1F *hRw1CentBJetEff = new TH1F("hRw1CentBJetEff","hRw1CentBJetEff",10,0,100);

  TH1F *hRw2CentJet = new TH1F("hRw2CentJet","hRw2CentJet",10,0,100);
  TH1F *hRw2CentBJet = new TH1F("hRw2CentBJet","hRw2CentBJet",10,0,100);
  TH1F *hRw2CentBJetEff = new TH1F("hRw2CentBJetEff","hRw2CentBJetEff",10,0,100);

  TH1F *hRw3CentJet = new TH1F("hRw3CentJet","hRw3CentJet",10,0,100);
  TH1F *hRw3CentBJet = new TH1F("hRw3CentBJet","hRw3CentBJet",10,0,100);
  TH1F *hRw3CentBJetEff = new TH1F("hRw3CentBJetEff","hRw3CentBJetEff",10,0,100);
  
  TH1F *hREtaJet = new TH1F("hREtaJet","hREtaJet",10,-etacut,etacut);
  TH1F *hREtaBJet = new TH1F("hREtaBJet","hREtaBJet",10,-etacut,etacut);
  TH1F *hREtaBJetEff = new TH1F("hREtaBJetEff","hREtaBJetEff",10,-etacut,etacut);
  
  TH1F *hRw1EtaJet = new TH1F("hRw1EtaJet","hRw1EtaJet",10,-etacut,etacut);
  TH1F *hRw1EtaBJet = new TH1F("hRw1EtaBJet","hRw1EtaBJet",10,-etacut,etacut);
  TH1F *hRw1EtaBJetEff = new TH1F("hRw1EtaBJetEff","hRw1EtaBJetEff",10,-etacut,etacut);
  
  TH1F *hRw2EtaJet = new TH1F("hRw2EtaJet","hRw2EtaJet",10,-etacut,etacut);
  TH1F *hRw2EtaBJet = new TH1F("hRw2EtaBJet","hRw2EtaBJet",10,-etacut,etacut);
  TH1F *hRw2EtaBJetEff = new TH1F("hRw2EtaBJetEff","hRw2EtaBJetEff",10,-etacut,etacut);
  
  TH1F *hRw3EtaJet = new TH1F("hRw3EtaJet","hRw3EtaJet",10,-etacut,etacut);
  TH1F *hRw3EtaBJet = new TH1F("hRw3EtaBJet","hRw3EtaBJet",10,-etacut,etacut);
  TH1F *hRw3EtaBJetEff = new TH1F("hRw3EtaBJetEff","hRw3EtaBJetEff",10,-etacut,etacut);
  
  
  TH1F *hRPtJet = new TH1F("hRPtJet","hRPtJet",10,100.,200.);
  TH1F *hRPtBJet = new TH1F("hRPtBJet","hRPtBJet",10,100.,200.);
  TH1F *hRPtBJetEff = new TH1F("hRPtBJetEff","hRPtBJetEff",10,100.,200.);

  TH1F *hRw1PtJet = new TH1F("hRw1PtJet","hRw1PtJet",10,100.,200.);
  TH1F *hRw1PtBJet = new TH1F("hRw1PtBJet","hRw1PtBJet",10,100.,200.);
  TH1F *hRw1PtBJetEff = new TH1F("hRw1PtBJetEff","hRw1PtBJetEff",10,100.,200.);

    TH1F *hRw2PtJet = new TH1F("hRw2PtJet","hRw2PtJet",10,100.,200.);
  TH1F *hRw2PtBJet = new TH1F("hRw2PtBJet","hRw2PtBJet",10,100.,200.);
  TH1F *hRw2PtBJetEff = new TH1F("hRw2PtBJetEff","hRw2PtBJetEff",10,100.,200.);

    TH1F *hRw3PtJet = new TH1F("hRw3PtJet","hRw3PtJet",10,100.,200.);
  TH1F *hRw3PtBJet = new TH1F("hRw3PtBJet","hRw3PtBJet",10,100.,200.);
  TH1F *hRw3PtBJetEff = new TH1F("hRw3PtBJetEff","hRw3PtBJetEff",10,100.,200.);

  
  
  hw1CentBJetEff->SetLineColor(2);
  hw1EtaBJetEff->SetLineColor(2);
  hw1PtBJetEff->SetLineColor(2);

  hw2CentBJetEff->SetLineColor(3);
  hw2EtaBJetEff->SetLineColor(3);
  hw2PtBJetEff->SetLineColor(3);

  hw3CentBJetEff->SetLineColor(4);
  hw3EtaBJetEff->SetLineColor(4);
  hw3PtBJetEff->SetLineColor(4);

  hw1CentBJetEff->SetMarkerColor(2);
  hw1EtaBJetEff->SetMarkerColor(2);
  hw1PtBJetEff->SetMarkerColor(2);

  hw2CentBJetEff->SetMarkerColor(3);
  hw2EtaBJetEff->SetMarkerColor(3);
  hw2PtBJetEff->SetMarkerColor(3);

  hw3CentBJetEff->SetMarkerColor(4);
  hw3EtaBJetEff->SetMarkerColor(4);
  hw3PtBJetEff->SetMarkerColor(4);

  hw1CentBJetEff->SetMarkerStyle(24);
  hw1EtaBJetEff->SetMarkerStyle(24);
  hw1PtBJetEff->SetMarkerStyle(24);

  hw2CentBJetEff->SetMarkerStyle(25);
  hw2EtaBJetEff->SetMarkerStyle(25);
  hw2PtBJetEff->SetMarkerStyle(25);

  hw3CentBJetEff->SetMarkerStyle(26);
  hw3EtaBJetEff->SetMarkerStyle(26);
  hw3PtBJetEff->SetMarkerStyle(26);


  hRw1CentBJetEff->SetLineColor(2);
  hRw1EtaBJetEff->SetLineColor(2);
  hRw1PtBJetEff->SetLineColor(2);
  
  hRw2CentBJetEff->SetLineColor(3);
  hRw2EtaBJetEff->SetLineColor(3);
  hRw2PtBJetEff->SetLineColor(3);

  hRw3CentBJetEff->SetLineColor(4);
  hRw3EtaBJetEff->SetLineColor(4);
  hRw3PtBJetEff->SetLineColor(4);

  hRw1CentBJetEff->SetMarkerColor(2);
  hRw1EtaBJetEff->SetMarkerColor(2);
  hRw1PtBJetEff->SetMarkerColor(2);

  hRw2CentBJetEff->SetMarkerColor(3);
  hRw2EtaBJetEff->SetMarkerColor(3);
  hRw2PtBJetEff->SetMarkerColor(3);

  hRw3CentBJetEff->SetMarkerColor(4);
  hRw3EtaBJetEff->SetMarkerColor(4);
  hRw3PtBJetEff->SetMarkerColor(4);

    hRw1CentBJetEff->SetMarkerStyle(24);
  hRw1EtaBJetEff->SetMarkerStyle(24);
  hRw1PtBJetEff->SetMarkerStyle(24);

  hRw2CentBJetEff->SetMarkerStyle(25);
  hRw2EtaBJetEff->SetMarkerStyle(25);
  hRw2PtBJetEff->SetMarkerStyle(25);

  hRw3CentBJetEff->SetMarkerStyle(26);
  hRw3EtaBJetEff->SetMarkerStyle(26);
  hRw3PtBJetEff->SetMarkerStyle(26);


  
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += nt->GetEntry(i);

      if(jtpt1<100.0) continue;
      if(discr_csvV1_1<0.9)continue;

      if(!hltCaloJet60 && !hltCaloJet80) continue;
      
      hCentJet->Fill(bin/2.);
      hEtaJet->Fill(jteta1);
      hPtJet->Fill(jtpt1);

      hRCentJet->Fill(bin/2.);
      hREtaJet->Fill(jteta1);
      hRPtJet->Fill(jtpt1);

      bool isTrig= false;

      if(hltCSV60&&hltCaloJet60) isTrig=true;
      else if(hltCSV80&&hltCaloJet80) isTrig=true;

      if(!isTrig) continue;
      
      hCentBJet->Fill(bin/2.);
      hEtaBJet->Fill(jteta1);
      hPtBJet->Fill(jtpt1);

      hRCentBJet->Fill(bin/2.);
      hREtaBJet->Fill(jteta1);
      hRPtBJet->Fill(jtpt1);
     
   }

   hCentBJetEff->Divide(hCentBJet,hCentJet,1.,1.,"B");
   hEtaBJetEff->Divide(hEtaBJet,hEtaJet,1.,1.,"B");
   hPtBJetEff->Divide(hPtBJet,hPtJet,1.,1.,"B");

   hRCentBJetEff->Divide(hRCentBJet,hRCentJet,1.,1.,"B");
   hREtaBJetEff->Divide(hREtaBJet,hREtaJet,1.,1.,"B");
   hRPtBJetEff->Divide(hRPtBJet,hRPtJet,1.,1.,"B");


   
   TF1 *fitCent = new TF1("fitCent","[0]*TMath::Erf((x-[1])/[2])",0.,200.);
   fitCent->SetParameters(0.5,20.5,30);
   hCentBJetEff->Fit(fitCent);



   double totCorr=0;
   double timesApplied=0;
  
   for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
     
     if(jtpt1<100.0) continue;
     if(discr_csvV1_1<0.9)continue;
     
     if(!hltCaloJet60 && !hltCaloJet80) continue;

      bool isTrig= false;

      if(hltCSV60&&hltCaloJet60) isTrig=true;
      else if(hltCSV80&&hltCaloJet80) isTrig=true;

      if(!isTrig) continue;

      double w1 = 1./fitCent->Eval(bin)/1.17068;
      double totW = w1;
      
     totCorr+=totW;
     timesApplied++;
      
      hw1CentBJet->Fill(bin/2.,totW);
      hw1EtaBJet->Fill(jteta1,totW);
      hw1PtBJet->Fill(jtpt1,totW);

      hRw1CentBJet->Fill(bin/2.,totW);
      hRw1EtaBJet->Fill(jteta1,totW);
      hRw1PtBJet->Fill(jtpt1,totW);
     
   }


   cout<<"average correction =  "<<totCorr/timesApplied<<endl;
   
   hw1CentBJetEff->Divide(hw1CentBJet,hCentJet,1.,1.,"B");
   hw1EtaBJetEff->Divide(hw1EtaBJet,hEtaJet,1.,1.,"B");
   hw1PtBJetEff->Divide(hw1PtBJet,hPtJet,1.,1.,"B");

   hRw1CentBJetEff->Divide(hRw1CentBJet,hRCentJet,1.,1.,"B");
   hRw1EtaBJetEff->Divide(hRw1EtaBJet,hREtaJet,1.,1.,"B");
   hRw1PtBJetEff->Divide(hRw1PtBJet,hRPtJet,1.,1.,"B");


   hw1EtaBJetEff->Fit("pol4");

   TF1 *fitEta = (TF1*)hw1EtaBJetEff->GetFunction("pol4");
   fitEta->SetName("fitEta");


  totCorr=0;
  timesApplied=0;
  
   
   for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
     
     if(jtpt1<100.0) continue;
     if(discr_csvV1_1<0.9)continue;
     
     if(!hltCaloJet60 && !hltCaloJet80) continue;

      bool isTrig= false;

      if(hltCSV60&&hltCaloJet60) isTrig=true;
      else if(hltCSV80&&hltCaloJet80) isTrig=true;

      if(!isTrig) continue;

      double w1 = 1./fitCent->Eval(bin)/1.17068;
      double w2 = 1./fitEta->Eval(jteta1)/1.20232;
      
      double totW = w1*w2;
      
      hw2CentBJet->Fill(bin/2.,totW);
      hw2EtaBJet->Fill(jteta1,totW);
      hw2PtBJet->Fill(jtpt1,totW);

      hRw2CentBJet->Fill(bin/2.,totW);
      hRw2EtaBJet->Fill(jteta1,totW);
      hRw2PtBJet->Fill(jtpt1,totW);

     totCorr+=totW;
     timesApplied++;

      
   }

      cout<<"average correction =  "<<totCorr/timesApplied<<endl;

   hw2CentBJetEff->Divide(hw2CentBJet,hCentJet,1.,1.,"B");
   hw2EtaBJetEff->Divide(hw2EtaBJet,hEtaJet,1.,1.,"B");
   hw2PtBJetEff->Divide(hw2PtBJet,hPtJet,1.,1.,"B");

   hRw2CentBJetEff->Divide(hRw2CentBJet,hRCentJet,1.,1.,"B");
   hRw2EtaBJetEff->Divide(hRw2EtaBJet,hREtaJet,1.,1.,"B");
   hRw2PtBJetEff->Divide(hRw2PtBJet,hRPtJet,1.,1.,"B");



   hw2PtBJetEff->Fit("pol3");


   totCorr=0;
   timesApplied=0;

   TF1 *fitPt = (TF1*)hw2PtBJetEff->GetFunction("pol3");
   fitPt->SetName("fitPt");
   for (Long64_t i=0; i<nentries;i++) {
     nt->GetEntry(i);
     
     if(jtpt1<100.0) continue;
     if(discr_csvV1_1<0.9)continue;
     
     if(!hltCaloJet60 && !hltCaloJet80) continue;

      
      bool isTrig= false;

      if(hltCSV60&&hltCaloJet60) isTrig=true;
      else if(hltCSV80&&hltCaloJet80) isTrig=true;

      if(!isTrig) continue;

      double w1 = 1./fitCent->Eval(bin)/1.17068;
      double w2 = 1./fitEta->Eval(jteta1)/1.20232;
      double w3 = 1./fitPt->Eval(jtpt1)/1.20232;
      if(jtpt1>200) w3=1./fitPt->Eval(200)/1.20232;
      
      double totW = w1*w2*w3;
      
      hw3CentBJet->Fill(bin/2.,totW);
      hw3EtaBJet->Fill(jteta1,totW);
      hw3PtBJet->Fill(jtpt1,totW);

      hRw3CentBJet->Fill(bin/2.,totW);
      hRw3EtaBJet->Fill(jteta1,totW);
      hRw3PtBJet->Fill(jtpt1,totW);

     totCorr+=totW;
     timesApplied++;

      
   }


   cout<<"average correction =  "<<totCorr/timesApplied<<endl;
   
   hw3CentBJetEff->Divide(hw3CentBJet,hCentJet,1.,1.,"B");
   hw3EtaBJetEff->Divide(hw3EtaBJet,hEtaJet,1.,1.,"B");
   hw3PtBJetEff->Divide(hw3PtBJet,hPtJet,1.,1.,"B");

   hRw3CentBJetEff->Divide(hRw3CentBJet,hRCentJet,1.,1.,"B");
   hRw3EtaBJetEff->Divide(hRw3EtaBJet,hREtaJet,1.,1.,"B");
   hRw3PtBJetEff->Divide(hRw3PtBJet,hRPtJet,1.,1.,"B");
   /*
   if(rebin){
     
     hCentBJetEff->Rebin(rebin);
     hCentBJetEff->Scale(1./(float)rebin);
     hw1CentBJetEff->Rebin(rebin);
     hw1CentBJetEff->Scale(1./(float)rebin);
     hw2CentBJetEff->Rebin(rebin);
     hw2CentBJetEff->Scale(1./(float)rebin);
     hw3CentBJetEff->Rebin(rebin);
     hw3CentBJetEff->Scale(1./(float)rebin);
   }
   */


   TLegend *leg = new TLegend(0.5,0.2,0.8,0.45);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hRCentBJetEff,"Uncorrected","p");
   leg->AddEntry(hRw1CentBJetEff,"Cent. corrected","p");
   leg->AddEntry(hRw2CentBJetEff,"Cent.-#eta corrected","p");
   leg->AddEntry(hRw3CentBJetEff,"Cent.-#eta-p_{T} corrected","p");
   
   TCanvas *c1=new TCanvas("c1","c1",600,600);

   hRCentBJetEff->SetMinimum(0.);
   hRCentBJetEff->SetMaximum(1.);   
   hRCentBJetEff->SetXTitle("centrality bin");
   hRCentBJetEff->SetYTitle("Online Tagging Efficiency");
   hRCentBJetEff->Draw();
   hRw1CentBJetEff->Draw("same");
   hRw2CentBJetEff->Draw("same");
   hRw3CentBJetEff->Draw("same");
   fitCent->Draw("same");
   leg->Draw();

   TCanvas *c2=new TCanvas("c2","c2",600,600);
   
   hREtaBJetEff->SetMinimum(0.);
   hREtaBJetEff->SetMaximum(1.);
   hREtaBJetEff->SetXTitle("#eta");
   hREtaBJetEff->SetYTitle("Online Tagging Efficiency");
   hREtaBJetEff->Draw();
   hRw1EtaBJetEff->Draw("same");
   hRw2EtaBJetEff->Draw("same");
   hRw3EtaBJetEff->Draw("same");
   fitEta->SetLineColor(2);
   fitEta->Draw("same");
   leg->Draw();

   TCanvas *c3=new TCanvas("c3","c3",600,600);
   
   hRPtBJetEff->SetMinimum(0.);
   hRPtBJetEff->SetMaximum(1.);
   hRPtBJetEff->SetXTitle("p_{T} (GeV)");
   hRPtBJetEff->SetYTitle("Online Tagging Efficiency");
   hRPtBJetEff->Draw();
   hRw1PtBJetEff->Draw("same");
   hRw2PtBJetEff->Draw("same");
   hRw3PtBJetEff->Draw("same");
   fitPt->SetLineColor(3);
   fitPt->Draw("same");
   leg->Draw();
   
   TFile *fout=NULL;
   if(write) {
     fout=new TFile("../correctionfiles/trigEffCorr.root","recreate");

     fitCent->SetParameters(fitCent->GetParameter(0), fitCent->GetParameter(1), fitCent->GetParameter(2));
     fitPt->Write();
     fitEta->Write();
     fitCent->Write();
     fout->Close();

     c1->SaveAs("onlineTrigEffCent.pdf");
     c1->SaveAs("onlineTrigEffCent.gif");
     c1->SaveAs("onlineTrigEffCent.C");
     c2->SaveAs("onlineTrigEffEta.pdf");
     c2->SaveAs("onlineTrigEffEta.gif");
     c2->SaveAs("onlineTrigEffEta.C");
     c3->SaveAs("onlineTrigEffPt.pdf");
     c3->SaveAs("onlineTrigEffPt.gif");
     c3->SaveAs("onlineTrigEffPt.C");

   }
   
}



