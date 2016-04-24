# ntuplizer

Builds compact ntuples from forest files.

## Order of operation

 1. Build data ntuples: ` sh buildtuples1.sh `
 2. Build qcd and b-jet filtered ntuples (uses data ntuples for centrality/vertex reweighting): ` sh buildtuples2.sh `
 3. Build and merge FCR ntuples (needs b-jet filtered ntuples): ` sh buildtuples3.sh `

Build everything altogether:

```
  sh run.sh &
```