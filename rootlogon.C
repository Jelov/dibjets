#include "../helpers/InitStyle.h"

void rootlogon()
{
  cout<<"root logon! "<<endl;
  InitStyle();
  
  gErrorIgnoreLevel = kWarning;
  gSystem->SetAclicMode(TSystem::kOpt);
  TH1::SetDefaultSumw2(true);
}
