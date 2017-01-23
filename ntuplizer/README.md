# ntuplizer

Builds compact ntuples from forest files.

## Order of operation

 1. Build data ntuples: ` sh buildtuples1.sh `
 2. Build qcd and b-jet filtered ntuples (uses data ntuples for centrality/vertex reweighting): ` sh buildtuples2.sh `
 3. Build and merge FCR ntuples (needs b-jet filtered ntuples): ` sh buildtuples3.sh `

Build everything altogether:

```bash
sh run.sh &
```


## Code description

Data processing is done in `buildtupledata.C`, Monte-Carlo in `buildtuplemc.C`. They are two very similar pieces of code which actually could have been merged…

The data is processed first (`buildtuples1.sh`), because the centrality and vertex distributions are built from data. Then QCD and b-jet filtered samples are processed (`buildtuples2.sh`). Finally, FCR samples are built and merged with b-jet filtered samples  (`buildtuples3.sh`) into `bfa`-files.

Both data and MC codes look for forest files under `samplesfolder`/`subfoldernames` and place the output files to `outputfolder` :

```C++
TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
TString samplesfolder = "/data_CMS/cms/mnguyen/bJet2015/mc";
```

Every ntuplizer task has the `code` which looks like: `mcPbqcdakPu4PF_djt.root`

- Letters 1-2: `mc`=MC, `dt`=Data
- Letters 3-4: `pp`=pp, `Pb`=PbPb (I should have used `HI` instead)
- Letters 5-7: sample code. It can be 
  - `qcd`=QCD, 
  - `bjt`=b-jet filtered sample
  - `bfc`=FCR samples
  - `bfa`=merged bjt+bfc (b-flavor augmented)
- rest of letters before underscore is the jet algorithm: `ak4PF`, `akPu4PF`, etc.
- After underscore - ntuple type: 
  - in `djt`  ntuple every entry = 1 event, and there are branches for jtpt1, jtpt2, etc. This is the file to draw dijet imbalance from.
  - in `inc` ntuple every jet is an entry.

If there is a need to add an additional sample, one has to create a code and add corresponding `if` statement to the `Init()` function. 

For example, MC PbPb sample with code `abc` and akPu4PF jet algorithm would be `mcPbabcakPu4PF` and corresponding if:

```C++
  if (PbPb && sample=="qbc") {
    samplesfolder+="PbPb/pythia6/";
    subfoldernames={"abc30","abc50","abc80","abc100","abc120"};
    pthats = {           30,     50,     80,     100,      120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  }
```

here the macro looks for `samplesfolder/PbPb/pythia6/abc{30,50,80,100,120}` folders and HiForest* files inside. Each folder corresponds to pT,hat(`pthats`) and cross-section (`CS`). Everything else should stay untouched. 

Run as: 

`root -l -q -b buildtuplemc.C+\(\"mcPbabcakPu4PF\"\) &`

Look for `mcPbabcakPu4PF_djt.root` and `mcPbabcakPu4PF_inc.root` in the `outputfolder`.



### Signal definition

Current definition assumes subid==0 and refpt>20. Jets which don't pass this requirement are considered combinatorial background.

```C++
auto isSignal = [&] (int N) {return subid[N]==0 && refpt[N]>20;};
```

### Variable definitions

Most of variables are self-explanatory. There are several groups of variables:

- `subid1:refpt1:rawpt1:jtpt1:...` connected to leading jet
- `subid2:refpt2:rawpt2:jtpt2:...`connected to subleading jet
- `subid3:refpt3:rawpt3:jtpt3:...`connected to third jet
- `SLord:subidSL:refptSL:rawptSL:jtptSL:…` skip-light jet, second *tagged* jet in the event, `SLord` is the number of the jet
- `NSLord:subidNSL:refptNSL:rawptNSL:jtptNSL` skip-light with negative tagger
- `SBord:subidSB:refptSB:rawptSB:jtptSB:…` second b-jet (defined by `abs(refparton_flavorForB)==5`)
- `Signal2ord:subidSignal2:refptSignal2:…` second signal jet (see definition of signal above).
- `SignalSLord:subidSignalSL:refptSignalSL:… ` skip-light signal jet.


Helper variable:`dphi21=acos(cos(jtphi[ind2]-jtphi[ind1]))`

`paircode` = integer value to represent type of the pair: 0 = BB, 1=CC, 2=BC/CB, 3=BX/XB, 4=XX/CX/XC, where X = not B or C, see `getPairCode()`.