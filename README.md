#Analysis procedure

##Remarks
* To prevent confusion, centrality bin number in the code is always 0-200 (except tagging efficiency correction functions), and shown as 0-100% only in the plotting.
* For additional help on ntuplizer see `ntuplizer/README.md`.
* For the additional help on the helpers and general code constructions (`Draw`, `Fill`, `macro m()`) , look at `helpers/README.md`.



## How to run the thing

* change working dir in helpers/config.h

```C++
  TString workdir = "/data_CMS/cms/lisniak/bjet2015/";
  TString tuplesfolderPbPb = workdir+"redoana/";
  TString tuplesfolderpp = workdir+"redoana/";
```

* create BXmistag functions:

```sh
mkdir correctionfiles
cd systematics
root -l -b -q BXmistags.C
```

* create trigger efficiency corrections (`iterTrigEffCorr.C`) and offline tagging corrections (`deriveOffTagEff.C`):

```sh
cd taggingcorrections
root -l -b -q  iterTrigEffCorr.C
sh runcorr.sh
cp -r tagging_jtsignal2v4_0p90/ ../correctionfiles
```

* tell me the truth!
  * `runall.sh` calculates all the xJ distributions with corrections applied for normal data (default) and smeared data (ppsmear1,2,3). You can look for distribution in `results_redoana_default/xJdphi.root` and mean values are stored in  `results_redoana_default/results.txt` .
  * `moneyplot.C` draws a few money plots
  * `plotXJ.C` draws beautiful xJ distributions. Make sure it reads the right input `results_***` folder!

```sh
cd asymmetry
sh runall.sh
#to check the progress
tail -f redoana_default
#after it finishes
root -l -b -q moneyplot.C\(\"redoana\"\)
root -l -b -q plotXJ.C
```

