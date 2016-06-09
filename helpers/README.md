# Helpers usage

## Plotting

The plotting routines are organised in the `helpers/plotting.h` header.

 - In the beginning of the macro **in the main macro function**

```
  macro m("macroname");
```

This will create a folder named `macroname` with all the subsequent plots and histograms.

This line creates object `m` which will be automatically desctructed in the end of the macro, and all the histograms will be saved into file at that moment. The constructor has a second boolean parameter `firstRun`, if this parameter is `false` then all histograms are read from this file instead of being filled with `Fill` function (see below), otherwise it's `true` by default.

```
  macro m("macroname",true); //default
  macro m("macroname",false); //histograms will be read from file
```
	

 - To create a histogram:

```
seth(nbins, xmin, xmax);
auto h1 = geth("histogramname1");
auto h2 = geth("histogramname2","histogramtitle;xaxistitle;yaxistitle");
```

`histogramtitle` is optional. After `seth` call all the histograms obtained with `geth` will have the same number of bins and x axis span. Also they will be saved in the root file in `macroname` folder. If `firstRun==false` then histogram is not created but read from the file.

- Normalize histogram

```
Normalize({h1,h2});
```
One can use as many histograms as needed inside `{}`. In order to normalize all the histograms created by `geth`, use `NormalizeAllHists();`.

- Draw histogram

In order to draw multiple histograms use

```
  Draw({h1,h2});
```
 
This function is designed to make quick and simple draw rather than sophisticated plot, but one can use many customizations to obtain more information, like:

```
plotputmean = true;
plotlegendpos = TopLeft;//TopRight,BottomRight,BottomLeft
plotymin = 0;
plotymax = 1;
...
```

In order to create a comparison plot with ratio/difference pad, use

```
DrawCompare(h1,h2);
```




## Loop ntuple
In order to simplify looping ntuples, the `Fill` command from `helpers/looptuple.h` is used across the analysis.
The syntax is the following:

```
Fill(file,[&] (dict &d) {
	//do all the work here      
  });
```

file has to be `TFile *` (for now) which contains tree (or ntuple) called `nt`. Inside the brackets `{}` one can write a code which will be executed for each entry of the tree. This allows to fill any histogram or perform any other type of analysis on the ntuple. 

In order to access values in the ntuple one has to use dictionary `d` as follows:

```
float v = d["jtpt"]; //reads value of jtpt branch
```

The full example is below. The following code loops through the ntuple, fills the histogram of xJ values, and plots it.

```
auto f = new TFile("dtppjpfak4PF_djt.root");

seth(10,0,1);
auto xj = geth("hxj");
Fill(f,[&] (dict &d) {

    if (d["jtpt1"]>100 && d["jtpt2"]>40 && d["dphi21"]>2.1)
      xj->Fill(d["jtpt2"]/d["jtpt1"],d["weight"]);
      
  });

Normalize({xj});
Draw({xj});

```