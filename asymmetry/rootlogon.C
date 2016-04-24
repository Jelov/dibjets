#include "../helpers/InitStyle.h"

void rootlogon()
{
  cout<<"root logon! "<<endl;
  InitStyle();
  
  TH1::SetDefaultSumw2(true);
}
