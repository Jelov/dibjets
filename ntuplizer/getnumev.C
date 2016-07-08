void getnumevpp()
{

  vector<TString> f = {
    "/data_CMS/cms/mnguyen/bJet2015/data//pp_PFLowPt/constSubV1_csvV2//merged_HiForestAOD.root",
"/data_CMS/cms/mnguyen/bJet2015/data//pp_PFHighPt/constSubV1_csvV2//merged_HiForestAOD.root"
  };

  int pp = 0;

  for (auto x:f) {
    cout<<x<<endl;

    auto file = new TFile(x);
    auto HiTree = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    HiTree->AddFriend("hltanalysis/HltTree");
    TString hlt = pp==0?"HLT_AK4PFJet60_Eta5p1_v1 && !HLT_AK4PFJet80_Eta5p1_v1":
      "HLT_AK4PFJet80_Eta5p1_v1";
    pp+=HiTree->GetEntries(hlt);

  }

  cout<<"pp "<<pp<<endl;
}

void getnumev()
{
  getnumevpp();

  vector<TString> f = {
    "/data_CMS/cms/mnguyen/bJet2015/data//PbPb_BJetSD/puTowerExclLimitV2/0000//merged_HiForestAOD.root",
    "/data_CMS/cms/mnguyen/bJet2015/data//PbPb_BJetSD/puTowerExclLimitV2/0001//merged_HiForestAOD_extra.root",
    "/data_CMS/cms/mnguyen/bJet2015/data//PbPb_BJetSD/puTowerExclLimitV2/0001//merged_HiForestAOD.root",
    "/data_CMS/cms/mnguyen/bJet2015/data//PbPb_BJetSD/puTowerExclLimitV2/0002//merged_HiForestAOD.root"};

  int central = 0;
  int midcent = 0;
  int periphl = 0;

  for (auto x:f) {
    cout<<x<<endl;

    auto file = new TFile(x);
    auto HiTree = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    central+=HiTree->GetEntries("hiBin<20");
    midcent+=HiTree->GetEntries("hiBin>=20 && hiBin<60");
    periphl+=HiTree->GetEntries("hiBin>=60");


  }

  cout<<"central "<<central<<endl;
  cout<<"midcent "<<midcent<<endl;
  cout<<"periphl "<<periphl<<endl;


}
