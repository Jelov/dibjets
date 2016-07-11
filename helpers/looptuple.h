#ifndef LOOPTUPLE_H
#define LOOPTUPLE_H


#include <functional>
#include <map>
#include "TNtuple.h"
#include "plotting.h"
#include "config.h"

//#define dict map<TString,float>

// list-absed dict - 66 s (use this one)
// my stupid map - 66 s
// std::map - 110 s...
class floatdict
{
public:
  vector<TString> keys;
  vector<float> values;
  void insert(TString key, float value)
  {
    auto p = std::find(keys.begin(), keys.end(), key);
    if (p == keys.end()) {
      keys.push_back(key);
      values.push_back(value);
    }
    else values[p-keys.begin()]=value;
  }
  //this method can insert in principle, but doesn't do that
  //to make sure values are not spoiled during reading
  bool brokenGet = false;
  float &operator[](const TString key)
  {
    auto p = std::find(keys.begin(), keys.end(), key);
    if (p == keys.end()) {
      cout<<"Key \""<<key<<"\" not found!"<<endl;
      brokenGet = true;
    }
    
    return values[p-keys.begin()];

  }

};

class dict
{
public:
  vector<TString> keys;
  vector<TTreeReaderValue<float> *> values;
  bool brokenGet = false;
  float &operator[](const TString key)
  {
    auto p = std::find(keys.begin(), keys.end(), key);
    if (p == keys.end()) {
      cout<<"Key \""<<key<<"\" not found!"<<endl;
      brokenGet = true;
    }
    
    return *(*values[p-keys.begin()]);

  }

};


class mystupidmap
{
  const static int maxN = 1<<16;
  TString keys[maxN];
  float values[maxN];
public:
  float &operator[](const TString key)
  {
    unsigned h = key.Hash() % maxN;//hashstring(key);
    keys[h] = key;
    return values[h];
  }
};

bool looptupledryrun = false;

void Fill(TFile *f, vector<TString> varsNeeded, const std::function<void(floatdict &)> & func, float processFraction = 1)
{
  if (looptupledryrun) processFraction = 0.01;
  if (!firstRunMacro) {
    cout<<" histograms have been read from file, skipping Fill function."<<endl;
    return;
  }

  //  for (auto x:varsNeeded) cout<<"\""<<x<<"\",";  cout<<endl;
	TTreeReader reader("nt",f);
	vector<TTreeReaderValue<float> *> values (varsNeeded.size());
	for (unsigned i=0;i<varsNeeded.size();i++)
		values[i] = new TTreeReaderValue<float>(reader,varsNeeded[i]);

	cout<<"Processing file "<<f->GetName()<<endl;
  bool fillportion = processFraction != 1;

	floatdict v;
	int nev = reader.GetEntries(true);
	int onep = nev/100;
	int evCounter = 0;
	TTimeStamp t0;

	while (reader.Next()) {
	  if (fillportion && evCounter>processFraction*nev) break;
		evCounter++;
		if (evCounter%onep==0) {
			std::cout << std::fixed; TTimeStamp t1; 
			cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
		}
float a = 0;

		for (unsigned i=0;i<varsNeeded.size();i++)
		  v.insert(varsNeeded[i], *(*(values[i])));//[(TString)varsNeeded[i]] = *(*(values[i]));
		func(v);
    //cout<<(int)v.brokenGet<<endl;
    if (v.brokenGet) {
      cout<<"You tried to access key which is not in the list of required branches"<<endl;
      break;
    }
	}
	cout<<endl;
}

void Fill(TFile *f, const std::function<void(dict &)> & func, float processFraction = 1)
{
  if (looptupledryrun) processFraction = 0.01;
  if (!firstRunMacro) {
    cout<<" histograms have been read from file, skipping Fill function."<<endl;
    return;
  }

  dict v;

  auto nt = (TTree *)f->Get("nt");
  TObjArray *brlist = nt->GetListOfBranches();
  int Nbranches = brlist->GetEntries();

  v.keys.resize(Nbranches);
  v.values.resize(Nbranches);

  for(int i = 0; i < Nbranches; ++i) { 
    TString brname = brlist->At(i)->GetName();
    v.keys[i] = brname;
  }

  TTreeReader reader("nt",f);

  for (int i=0;i<Nbranches;i++)
    v.values[i] = new TTreeReaderValue<float>(reader,v.keys[i]);



  cout<<"Processing file "<<f->GetName()<<endl;
  bool fillportion = processFraction != 1;


  int nev = reader.GetEntries(true);
  int onep = nev/100;
  int evCounter = 0;
  TTimeStamp t0;

  while (reader.Next()) {
    if (fillportion && evCounter>processFraction*nev) break;
    evCounter++;
    if (evCounter%onep==0) {
      std::cout << std::fixed; TTimeStamp t1; 
      cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
    }


    func(v);

    if (v.brokenGet) {
      cout<<"You tried to access key which is not present in the tree!"<<endl;
      break;
    }
  }
  cout<<endl;
}

//stub from old Fill
void Fill(TFile *f, vector<TString> varsNeeded,  const std::function<void(dict &)> & func, float processFraction = 1) 
{
  varsNeeded.clear();
  Fill(f,func,processFraction);
}


vector<TString> concat(vector<TString> a, vector<TString> b)
{
	auto v = vector<TString>();
	v.insert(v.end(),a.begin(), a.end());
	v.insert(v.end(),b.begin(), b.end());
	return v;
}
 

void WriteToFile(TString filename, vector<TString> varnames, vector<float> values)
{
  if (varnames.size()==0) return;
  TString s=varnames[0];
  for (unsigned i=1;i<varnames.size();i++) s+=":"+varnames[i];

  auto f = new TFile(filename, "recreate");
  auto nt = new TNtuple("nt","nt",s);
  nt->Fill(&values[0]);
  nt->Write();
  f->Close();

}

void WriteToFile(TString filename, map<TString, float> m)
{
  if (m.size()==0) return;
  vector<TString> varnames;
  vector<float> values;

  for (auto i:m) {
    varnames.push_back(i.first);
    values.push_back(i.second);
  }


  WriteToFile(filename, varnames, values);

}

map<TString,float> ReadFromFile(TString filename)
{
  vector<TString> varnames;
  map<TString, float> res;

  TFile *f = new TFile(filename);
  auto nt = (TTree *)f->Get("nt");
  auto l = nt->GetListOfBranches();
  for(int i = 0; i < l->GetEntries(); ++i)
    varnames.push_back(TString(l->At(i)->GetName()));

  TTreeReader reader("nt",f);
  vector<TTreeReaderValue<float> *> values (varnames.size());
  for (unsigned i=0;i<varnames.size();i++)
    values[i] = new TTreeReaderValue<float>(reader,varnames[i]);

  int nev = reader.GetEntries(true);

  if (nev==0) {
    cout<<"No entries in the file "<<filename<<endl;
    return res;
  }
  if (nev>1) {
    cout<<"There are more than 1 entry in the file "<<filename<<" ("<<nev<<"). Exiting."<<endl;
    return res;
  }

  reader.Next();
  for (unsigned i=0;i<varnames.size();i++)
    res[varnames[i]] =*(*(values[i]));
  
  return res;


}

#endif

