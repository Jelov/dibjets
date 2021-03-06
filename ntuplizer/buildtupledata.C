/* buildtupledata.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames
 * change jettree if you want to use different jet algo
 * collision = {PbPbBJet/PbPb/pp}
*/

#include "../helpers/parsecode.h"
#include <boost/functional/hash.hpp>
#include "Corrections.h"
#include <numeric>

TString jettree;
vector<TString> subfoldernames;

bool PbPb = false;

bool mockSL = false;
TString datatype = "";

TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
TString samplesfolder="/data_CMS/cms/mnguyen/bJet2015/data/";

int maxrepeat = 1;
const int NaN = -999;

const float hiHFcut = 5500;
const float etacut = 1.5;


bool applyjec = true;
Corrections jec;

bool applysmearing = false;
Smearing spp;


TTree *GetTree(TFile *f, TString treename)
{
  TTree *t = (TTree *)f->Get(treename);
  //PbPb bjet pthat120 has only unsubtracted jets!!!!! 
  //TODO: figure out
  //  if (t==0) t = (TTree *)f->Get("ak4PFJetAnalyzer/t");
  return t;
}

std::unordered_set<size_t> processedevents;
bool SeenEventAlready(unsigned int run,unsigned int lumi,unsigned long long event)
{
  vector<unsigned long long> v = {(unsigned long long)run,(unsigned long long)lumi,(unsigned long long)event};
  size_t hash = boost::hash_range(v.begin(),v.end());
  bool f = processedevents.find(hash)!=processedevents.end();
  if (!f) processedevents.insert(hash);
  // if (f) cout<<"seen already "<<run<<" _ "<<lumi<<" _ "<<event<<endl;
  return f;
}

vector<TString> list_files(const char *dirname, const char *exp=".*HiForestAOD.*\\.root")
{
  vector<TString> names;
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.Contains(TRegexp(exp))) {
        names.push_back(TString(dirname)+"/"+fname);
      }
    }
  }
  if (names.size()==0) return {dirname};

  return names;
}


int getEvents(TString folder, TString condition)
{
  auto files = list_files(Form("%s/%s",samplesfolder.Data(),folder.Data()));
  double x0=0;

  for (auto f:files) {
    TFile *f0 = TFile::Open(f);
    TTree *t0 = GetTree(f0,"hltanalysis/HltTree");//jettree);
    x0 += t0->GetEntries(condition);
    f0->Close();
  }

  return x0;
}

vector<double> weights;

void calculateWeights()
{
  cout<<"Calculating weights"<<endl;

  TString lowsample = subfoldernames[0];
  TString highsample = subfoldernames[1];

  int Njt80 = getEvents(highsample, "HLT_AK4PFJet80_Eta5p1_v1");
  int Njt80And60 = getEvents(highsample, "HLT_AK4PFJet80_Eta5p1_v1 && HLT_AK4PFJet60_Eta5p1_v1");
  
  weights = {(double)Njt80/(double)Njt80And60, 1};
}

vector<float> calculateWeightsBjet(TString filenamedj)
{
  cout<<"Calculating b-jet trigger weights"<<endl;

  TFile *f0 = new TFile(filenamedj);
  auto nt = (TTree *)f0->Get("nt");

  int Njt80 = nt->GetEntries("triggermatched && hltCSV80");
  int Njt80and60 = nt->GetEntries("triggermatched && hltCSV80 && hltCSV60");

  vector<float> w = {(float)Njt80/Njt80and60, 1};
  return w;
}

vector<float> calculateWeightsCaloJet(TString filenamedj)
{
  cout<<"Calculating calo-jet trigger weights"<<endl;

  TFile *f0 = new TFile(filenamedj);
  auto nt = (TTree *)f0->Get("nt");

  float njet100 = nt->GetEntries("hltCaloJet100");
  float overlapjet60 = nt->GetEntries("hltCaloJet60 && hltCaloJet100")/njet100;
  float overlapjet80 = nt->GetEntries("hltCaloJet80 && hltCaloJet100")/njet100;

  vector<float> w = {1/overlapjet60,1/overlapjet80,1};
  return w;
}



double getweight(TString sample, int trig60, int trig80)
{
  if (sample == subfoldernames[0] && trig60 && !trig80) return weights[0]; //pp_PFLowPt"
  if (sample == subfoldernames[1] && trig80) return 1; //pp_PFHighPt"

  return 0;
}

float matchingDistance(float jtphi1, float jteta1, float jtphi2, float jteta2)
{
  return (jteta1-jteta2)*(jteta1-jteta2)/0.6/0.6 + (jtphi1-jtphi2)*(jtphi1-jtphi2)/0.3/0.3;
}

bool matches(float dist)
{
  return dist<1.0;
} 

int triggeredLeadingJetCSV(float leadjtphi, float leadjteta, vector<Double_t> &trigpt, vector<Double_t> &trigphi, vector<Double_t> &trigeta)
{

  vector<int> btagged;
  vector<float> btaggedDist;
  for (unsigned i=0;i<trigpt.size();i++)
    for (unsigned j=i+1;j<trigpt.size();j++)
      if (trigpt[i]==trigpt[j]) btagged.push_back(i);

  for (auto ind:btagged) btaggedDist.push_back(matchingDistance(leadjtphi, leadjteta, trigphi[ind], trigeta[ind]));
  int bestmatchB = std::min_element(btaggedDist.begin(), btaggedDist.end()) - btaggedDist.begin();

  if (btaggedDist.size()>0 && matches(btaggedDist[bestmatchB])) return bestmatchB;

  return -1;
}

int triggeredLeadingJetCalo(float leadjtpt, float leadjtphi, float leadjteta, vector<Double_t> &trigpt, vector<Double_t> &trigphi, vector<Double_t> &trigeta)
{
  //returns closest deta-dphi jet

  // vector<float> caloDist;
  // for (unsigned i=0;i<trigpt.size();i++) caloDist.push_back(matchingDistance(leadjtphi, leadjteta, trigphi[i], trigeta[i]));
  // int bestmatchCalo = std::min_element(caloDist.begin(), caloDist.end()) - caloDist.begin();
  // if (caloDist.size()>0 && matches(caloDist[bestmatchCalo])) return bestmatchCalo;

  // return -1;

  vector<float> caloDist;
  for (unsigned i=0;i<trigpt.size();i++) {
    float d = matchingDistance(leadjtphi, leadjteta, trigphi[i], trigeta[i]);
    // caloDist.push_back(matches(d) ? abs(leadjtpt-trigpt[i]) : 9999);
    caloDist.push_back(matches(d) ? trigpt[i] : -1);
  }

  int bestmatchCalo = std::max_element(caloDist.begin(), caloDist.end()) - caloDist.begin();
  if (caloDist.size()>0) return bestmatchCalo;

  return -1;
}


void Init(TString sample)
{
  if (!PbPb && sample=="jpf") {
    subfoldernames = {"pp_PFLowPt/constSubV1_csvV2","pp_PFHighPt/constSubV1_csvV2"};
    calculateWeights();
  }
  else if (PbPb && sample=="bjt") {
    subfoldernames = {"PbPb_BJetSD/puTowerExclLimitV2/0000","PbPb_BJetSD/puTowerExclLimitV2/0001","PbPb_BJetSD/puTowerExclLimitV2/0002"};
  }
  else if (PbPb && sample=="j40") {
    subfoldernames = {"PbPb_Jet40/puTowerExclLimitV2/0000","PbPb_Jet40/puTowerExclLimitV2/0001","PbPb_Jet40/puTowerExclLimitV2/0002"};
    weights = {1.,1.,1.};
  }
  else if (PbPb && sample=="j60") {
    subfoldernames = {"PbPb_Jet6080/puTowerExclLimitV2/0000","PbPb_Jet6080/puTowerExclLimitV2/0001","PbPb_Jet6080/puTowerExclLimitV2/0002"};
    weights = {1.,1.,1.};
  }
  else if (PbPb && sample=="j80") {
    subfoldernames = {"PbPb_Jet6080/puTowerExclLimitV2/0000","PbPb_Jet6080/puTowerExclLimitV2/0001","PbPb_Jet6080/puTowerExclLimitV2/0002"};
    weights = {1.,1.,1.};
  }
  else if (PbPb && sample=="jcl") {
    subfoldernames = {"PbPb_Jet6080/puTowerExclLimitV2/0000","PbPb_Jet6080/puTowerExclLimitV2/0001","PbPb_Jet6080/puTowerExclLimitV2/0002",
                "PbPb_Jet100/puTowerExclLimitV2/0000","PbPb_Jet100/puTowerExclLimitV2/0001","PbPb_Jet100/puTowerExclLimitV2/0002"};
    weights = {1.,1.,1.,1.,1.,1.};
  }
  else cout<<"Don\'t know collision type: PbPb"<<PbPb<<", sample "<<sample<<endl;
}


void updatePbPbBtriggerweight(TString filename, vector<float> w)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float csv60, csv80;
  float triggermatched;
  float weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);

  nt->SetBranchAddress("hltCSV60",&csv60);
  nt->SetBranchAddress("hltCSV80",&csv80);
  nt->SetBranchAddress("triggermatched",&triggermatched);
  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    weight = 0;
    if (triggermatched && csv80) weight = w[1];
    if (triggermatched && csv60 && !csv80) weight = w[0];

    bw->Fill();
  }

  nt->Write("nt",TObject::kOverwrite);
  f->Close();


}

void updatePbPbCaloJetTriggerWeight(TString filename, vector<float> w)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float calojet60, calojet80,calojet100;
  float triggerPt;

  float weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);

  nt->SetBranchAddress("hltCaloJet60",&calojet60);
  nt->SetBranchAddress("hltCaloJet80",&calojet80);
  nt->SetBranchAddress("hltCaloJet100",&calojet100);
  nt->SetBranchAddress("triggerPt",&triggerPt);
  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    weight = 0;
    if (calojet100 && triggerPt>100)                weight = w[2];
    if (calojet80 && triggerPt>80 && triggerPt<=100) weight = w[1];
    if (calojet60 && triggerPt>60 && triggerPt<=80) weight = w[0];

    bw->Fill();
  }

  nt->Write("nt",TObject::kOverwrite);
  f->Close();


}


void updateweight(TString filename)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float prew, weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);
  nt->SetBranchAddress("prew",&prew);

  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);
    weight = prew;
    bw->Fill();
  }

  nt->Write("nt",TObject::kOverwrite);
  f->Close();


}

TF1 *fpp = 0, *fPb1, *fPb2, *fPb3;

void loadmockSLfunc()
{
  //not working for now b/c of the root version incompatibility
  // auto file = new TFile("BXmistagfunc.root");
  // fpp = (TF1 *)file->Get("fpp");
  // fPb1 = (TF1 *)file->Get("fPb1");
  // fPb2 = (TF1 *)file->Get("fPb2");
  // fPb3 = (TF1 *)file->Get("fPb3");

  fpp = new TF1("fpp","expo",40,200);
  fPb1 = new TF1("fPb1","expo",40,200);
  fPb2 = new TF1("fPb2","expo",40,200);
  fPb3 = new TF1("fPb3","expo",40,200);
  fpp->SetParameters(-4.94207,-0.00241127);
  fPb1->SetParameters(-4.01644,-0.0114637);
  fPb2->SetParameters(-4.6333,-0.00971588);
  fPb3->SetParameters(-5.03225,-0.00564332);


}

//condition on SL jet whether it's true SL (csv>0.9) or randomized using P(tag|L)
bool SLcondition(float csv, float pt, float bin)
{
  if (!mockSL) return csv>0.9;

  //mock SL!
  if (fpp == 0) loadmockSLfunc();

  //consider only csv<0.5 for mockSL
  if (csv>0.5) return false;

  float r = gRandom->Uniform();

  if (!PbPb) return r<fpp->Eval(pt);

  if (bin<20)            return r<fPb1->Eval(pt);
  if (bin>=20 && bin<60) return r<fPb2->Eval(pt);
  if (bin>=60)           return r<fPb3->Eval(pt);

  cout << "Should never get to this point" << endl;
  return false;
}

float getHighestTriggerPt(int CaloJet40, int CaloJet60, int CaloJet80, int CaloJet100, vector<double> &calo40pt, vector<double> &calo60pt,vector<double> &calo80pt, vector<double> &calo100pt)
{
  float triggerPt = NaN;
  if (CaloJet40)
    for (double x:calo40pt)
      if (x>triggerPt) triggerPt = x;

  if (CaloJet60)
    for (double x:calo60pt)
      if (x>triggerPt) triggerPt = x;

  if (CaloJet80)
    for (double x:calo80pt)
      if (x>triggerPt) triggerPt = x;

  if (CaloJet100)
    for (double x:calo100pt)
      if (x>triggerPt) triggerPt = x;

    return triggerPt;
}

float getcorrected(float pt, float eta, int bin)
{
  if (PbPb  && applyjec) return pt*(float)jec.factor(pt,eta,bin);
  if (!PbPb && applysmearing) return pt+(float)spp.rollpp(bin);

  return pt;
}

TF1 *fppbin = 0;
void makeppbin(TString ppcode)
{
  if (ppcode=="pp") {
    applysmearing = false;
    return;
  }

  applysmearing = true;
  float binmin = 0, binmax = 0;
  if (datatype=="p1") {binmin = 0; binmax = 20;}
  if (datatype=="p2") {binmin = 20; binmax = 60;}
  if (datatype=="p3") {binmin = 60; binmax = 200;}

  fppbin = new TF1("fppbin","[0]*exp(-[1]*x-[2]*x*x-[3]*x*x*x)",binmin,binmax);
  fppbin->SetParameters(4.53928e+04, 1.98556e-02, -2.61975e-05, 6.71250e-07);


}


int getbin(int b)
{
  if (PbPb) return b;
  
  //pp
  if (!applysmearing) return 1;

  return fppbin->GetRandom();
}

void swapi(int &a, int &b)
{
  int c = a;
  a=b;
  b=c;
}

void swapf(float &a, float &b)
{
  float c = a;
  a=b;
  b=c;
}

vector<int> ordered(vector<float> ptcor)
{
  vector<int> indices(ptcor.size());
  std::iota(begin(indices), end(indices), 0);

  std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return ptcor[a] > ptcor[b]; }
    );
  return indices;
}

void buildtupledata(TString code)//(TString collision = "PbPbBJet", TString jetalgo = "akVs4PFJetAnalyzer")
{
  if (!dt(code)) { cout<<"Not data: "<<code<<", exiting..."<<endl; return;}

  Corrections jec;
  
  PbPb = isPbPb(code);
  TString sample = getSample(code);
  jettree = getjettree(code);
  mockSL = IsMockSL(code);
  datatype = getdatatype(code);

  makeppbin(datatype);

  if (mockSL) maxrepeat = 10;

  Init(sample);

  TString outputfilenamedj = outputfolder+"/"+code+"_djt.root";
  TString outputfilenameinc = outputfolder+"/"+code+"_inc.root";

  TString djvars = TString("run:lumi:event:prew:triggermatched:bin:vz:hiHF:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltCaloJet100:triggerPt:hltPFJet60:hltPFJet80:dijet:")+
      "hltCalo40jtpt:hltCalo40jtphi:hltCalo40jteta:hltCalo60jtpt:hltCalo60jtphi:hltCalo60jteta:hltCalo80jtpt:hltCalo80jtphi:hltCalo80jteta:hltCSV60jtpt:hltCSV60jtphi:hltCSV60jteta:hltCSV80jtpt:hltCSV80jtphi:hltCSV80jteta:"+
      "numTagged:"
      "rawpt1:jtpt1:jtptsansjec1:jtphi1:jteta1:discr_csvV1_1:ndiscr_csvV1_1:discr_ssvHighEff1:discr_ssvHighPur1:svtxm1:discr_prob1:svtxdls1:svtxpt1:svtxntrk1:nsvtx1:nselIPtrk1:"+
      "rawpt2:jtpt2:jtptsansjec2:jtphi2:jteta2:discr_csvV1_2:ndiscr_csvV1_2:discr_ssvHighEff2:discr_ssvHighPur2:svtxm2:discr_prob2:svtxdls2:svtxpt2:svtxntrk2:nsvtx2:nselIPtrk2:dphi21:"+
      "rawpt3:jtpt3:jtptsansjec3:jtphi3:jteta3:discr_csvV1_3:ndiscr_csvV1_3:discr_ssvHighEff3:discr_ssvHighPur3:svtxm3:discr_prob3:svtxdls3:svtxpt3:svtxntrk3:nsvtx3:nselIPtrk3:dphi31:dphi32:"+
      "SLord:rawptSL:jtptSL:jtptsansjecSL:jtphiSL:jtetaSL:discr_csvV1_SL:ndiscr_csvV1_SL:discr_ssvHighEffSL:discr_ssvHighPurSL:svtxmSL:discr_probSL:svtxdlsSL:svtxptSL:svtxntrkSL:nsvtxSL:nselIPtrkSL:dphiSL1:"+
      "NSLord:rawptNSL:jtptNSL:jtptsansjecNSL:jtphiNSL:jtetaNSL:discr_csvV1_NSL:ndiscr_csvV1_NSL:discr_ssvHighEffNSL:discr_ssvHighPurNSL:svtxmNSL:discr_probNSL:svtxdlsNSL:svtxptNSL:svtxntrkNSL:nsvtxNSL:nselIPtrkNSL:dphiNSL1";


  TString incvars = TString("prew:triggermatched:bin:vz:hiHF:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltCaloJet100:triggerPt:hltPFJet60:hltPFJet80:")+
      "rawpt:jtpt:jtptsansjec:jtphi:jteta:discr_csvV1:ndiscr_csvV1:discr_ssvHighEff:discr_ssvHighPur:svtxm:discr_prob:svtxdls:svtxpt:svtxntrk:nsvtx:nselIPtrk";


  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","ntdj",djvars);

  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","ntinc",incvars);

  for (unsigned i=0;i<subfoldernames.size();i++) {
    //get all files for unmerged forests
    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = new TFile(filename);
    TString treename = jettree;//f->Get(jettree) != 0 ? jettree : "ak3PFJetAnalyzer";
    TTreeReader reader(treename,f);
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> rawpt(reader, "rawpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<float> discr_csvV1(reader, "discr_csvV1");
    TTreeReaderArray<float> ndiscr_csvV1(reader, "ndiscr_csvV1");

    TTreeReaderArray<float> discr_ssvHighEff(reader, "discr_ssvHighEff");
    TTreeReaderArray<float> discr_ssvHighPur(reader, "discr_ssvHighPur");

    TTreeReaderArray<float> discr_prob(reader, "discr_prob");
    TTreeReaderArray<float> svtxm(reader, "svtxm");
    TTreeReaderArray<float> svtxdls(reader, "svtxdls");
    TTreeReaderArray<float> svtxpt(reader, "svtxpt");

    TTreeReaderArray<int> svtxntrk(reader, "svtxntrk");
    TTreeReaderArray<int> nsvtx(reader, "nsvtx");
    TTreeReaderArray<int> nselIPtrk(reader, "nselIPtrk");

    TTreeReaderArray<float> *muMax=0, *muMaxTRK=0, *muMaxGBL=0;
    if (PbPb) {
      muMax = new TTreeReaderArray<float> (reader, "muMax");
      muMaxTRK = new TTreeReaderArray<float>(reader, "muMaxTRK");
      muMaxGBL = new TTreeReaderArray<float>(reader, "muMaxGBL");
    }


    //HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1 HLT_HIPuAK4CaloJet80_Eta5p1_v1

    TString calojet40trigger = !PbPb ? "HLT_AK4CaloJet40_Eta5p1_v1" : "HLT_HIPuAK4CaloJet40_Eta5p1_v1";
    TString calojet40triggerv2 = !PbPb ? "HLT_AK4CaloJet40_Eta5p1_v1" : "HLT_HIPuAK4CaloJet40_Eta5p1_v2";
    TString calojet60trigger = !PbPb ? "HLT_AK4CaloJet60_Eta5p1_v1" : "HLT_HIPuAK4CaloJet60_Eta5p1_v1";
    TString calojet80trigger = !PbPb ? "HLT_AK4CaloJet80_Eta5p1_v1" : "HLT_HIPuAK4CaloJet80_Eta5p1_v1";
    TString calojet100trigger = !PbPb ? "HLT_AK4CaloJet100_Eta5p1_v1" : "HLT_HIPuAK4CaloJet100_Eta5p1_v1";
    //dummy vars in PbPb case
    TString pfjet60trigger = !PbPb ? "HLT_AK4PFJet60_Eta5p1_v1" : "LumiBlock";
    TString pfjet80trigger = !PbPb ? "HLT_AK4PFJet80_Eta5p1_v1" : "LumiBlock";
    TString csv60trigger = !PbPb ? "HLT_AK4PFBJetBCSV60_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1";
    TString csv80trigger = !PbPb ? "HLT_AK4PFBJetBCSV80_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1";

    //PbPb pprimaryVertexFilter && pclusterCompatibilityFilter do nothing
    vector<TString> filterNames;
    if (PbPb) filterNames = {"pcollisionEventSelection", "HBHENoiseFilterResultRun2Loose"};
    else filterNames = {"pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "pBeamScrapingFilter"}; 

    TTreeReader readerhlt("hltanalysis/HltTree",f);
    TTreeReaderValue<int> PFJet60(readerhlt, pfjet60trigger);
    TTreeReaderValue<int> PFJet80(readerhlt, pfjet80trigger);


    TTreeReaderValue<int> CaloJet40(readerhlt, calojet40trigger);
    TTreeReaderValue<int> CaloJet40v2(readerhlt, calojet40triggerv2);
    TTreeReaderValue<int> CaloJet60(readerhlt, calojet60trigger);
    TTreeReaderValue<int> CaloJet80(readerhlt, calojet80trigger);
    TTreeReaderValue<int> CaloJet100(readerhlt, calojet100trigger);    

    TTreeReaderValue<int> CSV60(readerhlt, csv60trigger);
    TTreeReaderValue<int> CSV80(readerhlt, csv80trigger);

    TTreeReader readercsv60object("hltobject/HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v",f);
    TTreeReaderValue<vector<Double_t> > csv60pt(readercsv60object, "pt");
    TTreeReaderValue<vector<Double_t> > csv60eta(readercsv60object, "eta");
    TTreeReaderValue<vector<Double_t> > csv60phi(readercsv60object, "phi");
    
    TTreeReader readercsv80object("hltobject/HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v",f);
    TTreeReaderValue<vector<Double_t> > csv80pt(readercsv80object, "pt");
    TTreeReaderValue<vector<Double_t> > csv80eta(readercsv80object, "eta");
    TTreeReaderValue<vector<Double_t> > csv80phi(readercsv80object, "phi");
    
    TTreeReader readerCalo40object("hltobject/HLT_HIPuAK4CaloJet40_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo40pt(readerCalo40object, "pt");
    TTreeReaderValue<vector<Double_t> > calo40eta(readerCalo40object, "eta");
    TTreeReaderValue<vector<Double_t> > calo40phi(readerCalo40object, "phi");

    TTreeReader readerCalo60object("hltobject/HLT_HIPuAK4CaloJet60_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo60pt(readerCalo60object, "pt");
    TTreeReaderValue<vector<Double_t> > calo60eta(readerCalo60object, "eta");
    TTreeReaderValue<vector<Double_t> > calo60phi(readerCalo60object, "phi");
    
    TTreeReader readerCalo80object("hltobject/HLT_HIPuAK4CaloJet80_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo80pt(readerCalo80object, "pt");
    TTreeReaderValue<vector<Double_t> > calo80eta(readerCalo80object, "eta");
    TTreeReaderValue<vector<Double_t> > calo80phi(readerCalo80object, "phi");

    TTreeReader readerCalo100object("hltobject/HLT_HIPuAK4CaloJet100_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo100pt(readerCalo100object, "pt");
    TTreeReaderValue<vector<Double_t> > calo100eta(readerCalo100object, "eta");
    TTreeReaderValue<vector<Double_t> > calo100phi(readerCalo100object, "phi");

    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    TTreeReaderValue<int> bin(readerevt, "hiBin");
    TTreeReaderValue<float> hiHF(readerevt, "hiHF");
    
    TTreeReaderValue<unsigned int> run(readerevt, "run");
    TTreeReaderValue<unsigned int> lumi(readerevt, "lumi");
    TTreeReaderValue<unsigned long long> event(readerevt, "evt");

    TTreeReader readerskim("skimanalysis/HltTree",f);

    vector<TTreeReaderValue<int> *>filters;
    for (auto f:filterNames)
      filters.push_back(new TTreeReaderValue<int>(readerskim, f));
      
    cout<<"added filters"<<endl;
    
    int nev = reader.GetEntries(true); cout<<nev<<endl;
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;
    
    //allows repeating event maxrepeat times
    int repeatcounter = maxrepeat;
    bool continuereading = true;
    while (true) {
      repeatcounter++;
      if (repeatcounter>=maxrepeat) {
        continuereading = reader.Next();
        readerhlt.Next();
        readerevt.Next();
        readerskim.Next();
        readercsv60object.Next();
        readercsv80object.Next();
        readerCalo40object.Next();
        readerCalo60object.Next();
        readerCalo80object.Next();
        readerCalo100object.Next();
        repeatcounter = 0;
      }
      if (!continuereading) break;

      // for testing - only 2% of data
      // if (evCounter>2*onep) break;

      evCounter++;
      if (evCounter%onep==0) {
        std::cout << std::fixed;
        TTimeStamp t1; 
        cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      if (PbPb && SeenEventAlready(*run,*lumi,*event)) continue; //in pp weight is 0 if repeating

      int bPFJet60 = !PbPb ? *PFJet60 : 1;
      int bPFJet80 = !PbPb ? *PFJet80 : 1;

      //int jet40 = *CaloJet40 || *CaloJet40v2;

      float weight = 1;

      if (!PbPb)
        weight = getweight(subfoldernames[i], bPFJet60, bPFJet80);

      if (PbPb && sample=="j60")
        weight = *CaloJet60;

      if (PbPb && sample=="j80")
        weight = *CaloJet80;

      //good event is vertex cut and noise cuts
      bool goodevent = abs(*vz)<15 && *hiHF<hiHFcut;
      for (auto f:filters) 
        goodevent&=*(*f);

      if (!goodevent) weight = 0;

      if (weight==0) continue;


      int ind1=-1, ind2=-1, ind3=-1, indSL=-1, indNSL=-1; //indices of leading/subleading jets in jet array
      int indTrigCSV60=-1, indTrigCSV80=-1, indTrigCalo40=-1, indTrigCalo60=-1, indTrigCalo80=-1;
      int SLord = 0, NSLord = 0;
      bool foundJ1=false, foundJ2 = false, foundJ3 = false, foundSL = false, foundNSL = false; //found/not found yet, for convenience

      bool triggermatched = false; float triggerPt = NaN;
      int numTagged = 0;

      int cbin = getbin(*bin);

      vector<float> ptcor(*nref);
      for (int j=0;j<*nref;j++)
        ptcor[j] = getcorrected(jtpt[j],jteta[j],cbin);

      vector<int> orderedind = ordered(ptcor);
      
      if (goodevent)
        for (int k=0;k<*nref;k++) {
          int j=orderedind[k];
          //acceptance selection
          if (abs(jteta[j])>etacut) continue;
          //muon cuts
          if (PbPb) {
            if((*muMax)[j]/rawpt[j]>0.95) continue;
            if( ((*muMaxTRK)[j]-(*muMaxGBL)[j]) / ((*muMaxTRK)[j]+(*muMaxGBL)[j]) > 0.1) continue;
          }
  
          float jtptcorrected = ptcor[j];//getcorrected(jtpt[j],jteta[j],cbin);

          if (!foundJ1) { //looking for the leading jet
              ind1 = j;
              foundJ1=true;

	        if (PbPb) {
		         indTrigCSV60 = triggeredLeadingJetCSV(jtphi[j], jteta[j], *csv60pt, *csv60phi, *csv60eta);
		         indTrigCSV80 = triggeredLeadingJetCSV(jtphi[j], jteta[j], *csv80pt, *csv80phi, *csv80eta);
             // actually, I don't need matching of calo jet to leading jet
             indTrigCalo40 = triggeredLeadingJetCalo(ptcor[j],jtphi[j], jteta[j], *calo40pt, *calo40phi, *calo40eta);//jtpt[j]
             indTrigCalo60 = triggeredLeadingJetCalo(ptcor[j],jtphi[j], jteta[j], *calo60pt, *calo60phi, *calo60eta);//jtpt[j]
		         indTrigCalo80 = triggeredLeadingJetCalo(ptcor[j],jtphi[j], jteta[j], *calo80pt, *calo80phi, *calo80eta);//jtpt[j]

             triggerPt = getHighestTriggerPt(*CaloJet40,*CaloJet60, *CaloJet80, *CaloJet100, *calo40pt,*calo60pt,*calo80pt,*calo100pt);
	        }
             
	        triggermatched = !PbPb || indTrigCSV60!=-1 || indTrigCSV80!=-1;
	       } else
            if (foundJ1 && !foundJ2) {
              ind2 = j;
              foundJ2 = true;
            } else
            if (foundJ1 && foundJ2 && !foundJ3) {
              ind3 = j;
              foundJ3 = true;
            }

          //we need ordinal number of SL jets, so counting until found
          //indSL != SLord because some jets are not in acceptance region
            if (!foundSL) SLord++;

          //ind1!=j otherwise SL will be = J1
            if (foundJ1 && ind1!=j && !foundSL && SLcondition(discr_csvV1[j], ptcor[j], cbin)) {
              indSL = j;
              foundSL = true;
            }  

            if (!foundNSL) NSLord++;

            //ind1!=j otherwise SL will be = J1
            if (foundJ1 && ind1!=j && !foundNSL && ndiscr_csvV1[j]>0.9) {
              indNSL = j;
              foundNSL = true;
            }

            if (SLcondition(discr_csvV1[j], ptcor[j], cbin)) numTagged++;


            //at this point foundLJ = true always, so triggermatched is determined
            vector<float> vinc = {weight, (float)triggermatched, (float) cbin, *vz, *hiHF,(float)*CSV60, (float)*CSV80,
              (float)*CaloJet40, (float)*CaloJet60, (float)*CaloJet80,(float)*CaloJet100,triggerPt,
				  (float)bPFJet60,(float)bPFJet80, rawpt[j], jtptcorrected, jtpt[j], jtphi[j], jteta[j], discr_csvV1[j],ndiscr_csvV1[j],discr_ssvHighEff[j],discr_ssvHighPur[j],svtxm[j],discr_prob[j],
              svtxdls[j],svtxpt[j],(float)svtxntrk[j],(float)nsvtx[j],(float)nselIPtrk[j]};
    
            if (!mockSL) //no need for inc ntuple in mockSL mode
              ntinc->Fill(&vinc[0]);
        }

        

      // float jtpt1cor = ptcor[ind1];//getcorrected(jtpt[ind1],jteta[ind1],cbin);
      // float jtpt2cor = ptcor[ind2];//getcorrected(jtpt[ind2],jteta[ind2],cbin);

      // if (jtpt2cor>jtpt1cor) {swapi(ind1,ind2); swapf(jtpt1cor,jtpt2cor); cout<<"WHAAAT???"<<endl;}
	//{int c=ind2; ind2=ind1; ind1=c;}// it doesn't mean that jtpt3 cannot enter the game!

      //fill dijet ntuple
      vector<float> vdj;

      vdj = {(float)*run, (float)*lumi, (float)*event, weight, (float)triggermatched, (float)cbin, *vz,*hiHF,
        (float)*CSV60, (float)*CSV80,(float)*CaloJet40,(float)*CaloJet60, (float)*CaloJet80,(float)*CaloJet100,triggerPt,(float)bPFJet60,(float)bPFJet80, 
        foundJ1 && foundJ2 ? (float)1 : (float)0,

        indTrigCalo40!=-1 ? (float)(*calo40pt)[indTrigCalo40] : NaN,
        indTrigCalo40!=-1 ? (float)(*calo40phi)[indTrigCalo40] : NaN,
        indTrigCalo40!=-1 ? (float)(*calo40eta)[indTrigCalo40] : NaN,

        indTrigCalo60!=-1 ? (float)(*calo60pt)[indTrigCalo60] : NaN,
        indTrigCalo60!=-1 ? (float)(*calo60phi)[indTrigCalo60] : NaN,
        indTrigCalo60!=-1 ? (float)(*calo60eta)[indTrigCalo60] : NaN,

        indTrigCalo80!=-1 ? (float)(*calo80pt)[indTrigCalo80] : NaN,
        indTrigCalo80!=-1 ? (float)(*calo80phi)[indTrigCalo80] : NaN,
        indTrigCalo80!=-1 ? (float)(*calo80eta)[indTrigCalo80] : NaN,

        indTrigCSV60!=-1  ? (float)(*csv60pt)[indTrigCSV60] : NaN,
        indTrigCSV60!=-1  ? (float)(*csv60phi)[indTrigCSV60] : NaN,
        indTrigCSV60!=-1  ? (float)(*csv60eta)[indTrigCSV60] : NaN,

        indTrigCSV80!=-1  ? (float)(*csv80pt)[indTrigCSV80] : NaN,
        indTrigCSV80!=-1  ? (float)(*csv80phi)[indTrigCSV80] : NaN,
        indTrigCSV80!=-1  ? (float)(*csv80eta)[indTrigCSV80] : NaN,
                 
	     (float)numTagged,
                       
        foundJ1 ? rawpt[ind1] : NaN,
	      foundJ1 ? ptcor[ind1] : NaN, //getcorrected(jtpt[ind1],jteta[ind1],cbin) : NaN,
        foundJ1 ? jtpt[ind1] : NaN,
        foundJ1 ? jtphi[ind1] : NaN,
        foundJ1 ? jteta[ind1] : NaN,
        foundJ1 ? discr_csvV1[ind1] : NaN,
        foundJ1 ? ndiscr_csvV1[ind1] : NaN,
        foundJ1 ? discr_ssvHighEff[ind1] : NaN,
        foundJ1 ? discr_ssvHighPur[ind1] : NaN,
        foundJ1 ? svtxm[ind1] : NaN,
        foundJ1 ? discr_prob[ind1] : NaN,
        foundJ1 ? svtxdls[ind1] : NaN,
        foundJ1 ? svtxpt[ind1] : NaN,
        foundJ1 ? (float)svtxntrk[ind1] : NaN,
        foundJ1 ? (float)nsvtx[ind1] : NaN,
        foundJ1 ? (float)nselIPtrk[ind1] : NaN,

        foundJ2 ? rawpt[ind2] : NaN,
	      foundJ2 ? ptcor[ind2] : NaN, //getcorrected(jtpt[ind2],jteta[ind2],cbin) : NaN,
        foundJ2 ? jtpt[ind2] : NaN,
        foundJ2 ? jtphi[ind2] : NaN,
        foundJ2 ? jteta[ind2] : NaN,
        foundJ2 ? discr_csvV1[ind2] : NaN,
        foundJ2 ? ndiscr_csvV1[ind2] : NaN,
        foundJ2 ? discr_ssvHighEff[ind2] : NaN,
        foundJ2 ? discr_ssvHighPur[ind2] : NaN,
        foundJ2 ? svtxm[ind2] : NaN,
        foundJ2 ? discr_prob[ind2] : NaN,
        foundJ2 ? svtxdls[ind2] : NaN, 
        foundJ2 ? svtxpt[ind2] : NaN,
        foundJ2 ? (float)svtxntrk[ind2] : NaN,
        foundJ2 ? (float)nsvtx[ind2] : NaN,
        foundJ2 ? (float)nselIPtrk[ind2] : NaN,
        foundJ2 && foundJ1 ? acos(cos(jtphi[ind2]-jtphi[ind1])) : NaN,
    
        foundJ3 ? rawpt[ind3] : NaN,
        foundJ3 ? ptcor[ind3] : NaN, //getcorrected(jtpt[ind3],jteta[ind3],cbin) : NaN,
        foundJ3 ? jtpt[ind3] : NaN,
        foundJ3 ? jtphi[ind3] : NaN,
        foundJ3 ? jteta[ind3] : NaN,
        foundJ3 ? discr_csvV1[ind3] : NaN,
        foundJ3 ? ndiscr_csvV1[ind3] : NaN,
        foundJ3 ? discr_ssvHighEff[ind3] : NaN,
        foundJ3 ? discr_ssvHighPur[ind3] : NaN,
        foundJ3 ? svtxm[ind3] : NaN,
        foundJ3 ? discr_prob[ind3] : NaN,
        foundJ3 ? svtxdls[ind3] : NaN, 
        foundJ3 ? svtxpt[ind3] : NaN,
        foundJ3 ? (float)svtxntrk[ind3] : NaN,
        foundJ3 ? (float)nsvtx[ind3] : NaN,
        foundJ3 ? (float)nselIPtrk[ind3] : NaN,
        foundJ3 && foundJ1 ? acos(cos(jtphi[ind3]-jtphi[ind1])) : NaN,
        foundJ3 && foundJ2 ? acos(cos(jtphi[ind3]-jtphi[ind2])) : NaN,

        foundSL ? (float)SLord : NaN,
        foundSL ? rawpt[indSL] : NaN,
        foundSL ? ptcor[indSL] : NaN,//getcorrected(jtpt[indSL],jteta[indSL],cbin) : NaN,
        foundSL ? jtpt[indSL] : NaN,
        foundSL ? jtphi[indSL] : NaN,
        foundSL ? jteta[indSL] : NaN,
        foundSL ? discr_csvV1[indSL] : NaN,
        foundSL ? ndiscr_csvV1[indSL] : NaN,
        foundSL ? discr_ssvHighEff[indSL] : NaN,
        foundSL ? discr_ssvHighPur[indSL] : NaN,
        foundSL ? svtxm[indSL] : NaN,
        foundSL ? discr_prob[indSL] : NaN,
        foundSL ? svtxdls[indSL] : NaN, 
        foundSL ? svtxpt[indSL] : NaN,
        foundSL ? (float)svtxntrk[indSL] : NaN,
        foundSL ? (float)nsvtx[indSL] : NaN,
        foundSL ? (float)nselIPtrk[indSL] : NaN,
        foundSL && foundJ1 ? acos(cos(jtphi[indSL]-jtphi[ind1])) : NaN,

        foundNSL ? (float)NSLord : NaN,
        foundNSL ? rawpt[indNSL] : NaN,
        foundNSL ? ptcor[indNSL] : NaN,//getcorrected(jtpt[indNSL],jteta[indNSL],cbin) : NaN,
        foundNSL ? jtpt[indNSL] : NaN,
        foundNSL ? jtphi[indNSL] : NaN,
        foundNSL ? jteta[indNSL] : NaN,
        foundNSL ? discr_csvV1[indNSL] : NaN,
        foundNSL ? ndiscr_csvV1[indNSL] : NaN,
        foundNSL ? discr_ssvHighEff[indNSL] : NaN,
        foundNSL ? discr_ssvHighPur[indNSL] : NaN,
        foundNSL ? svtxm[indNSL] : NaN,
        foundNSL ? discr_prob[indNSL] : NaN,
        foundNSL ? svtxdls[indNSL] : NaN, 
        foundNSL ? svtxpt[indNSL] : NaN,
        foundNSL ? (float)svtxntrk[indNSL] : NaN,
        foundNSL ? (float)nsvtx[indNSL] : NaN,
        foundNSL ? (float)nselIPtrk[indNSL] : NaN,
        foundNSL && foundJ1 ? acos(cos(jtphi[indNSL]-jtphi[ind1])) : NaN


      };

      //==(mockSL && foundSL) || !mockSL
      if (!mockSL || foundSL)
        ntdj->Fill(&vdj[0]);



    }

    f->Close();
    }
  }

  foutdj->cd();
  ntdj->Write("nt",TObject::kOverwrite);
  foutdj->Close();

  foutinc->cd();
  ntinc->Write("nt",TObject::kOverwrite);
  foutinc->Close();

  cout<<endl;
  cout<<"Total input entries "<<totentries<<endl;

  //making centrality-dependent ntuples
  //PutInCbins(outputfolder, code, {{0,40}, {80,200}});

  if (PbPb && sample=="bjt"){
    auto w = calculateWeightsBjet(outputfilenamedj);

    updatePbPbBtriggerweight(outputfilenamedj,w);
    updatePbPbBtriggerweight(outputfilenameinc,w);
  } else if (PbPb && sample=="jcl"){
    auto w = calculateWeightsCaloJet(outputfilenamedj);

    updatePbPbCaloJetTriggerWeight(outputfilenamedj,w);
    updatePbPbCaloJetTriggerWeight(outputfilenameinc,w);
  }
  else {
    updateweight(outputfilenamedj);
    updateweight(outputfilenameinc);
  }

}
