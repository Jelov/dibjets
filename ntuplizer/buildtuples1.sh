root -l <<< ".L buildtupledata.C+"

#Build all data tuples
root -l -q -b buildtupledata.C+\(\"dtppjpfak4PF\"\) >out/dtppjpfak4PF &
root -l -q -b buildtupledata.C+\(\"dtppjpfak3PF\"\) >out/dtppjpfak3PF &

root -l -q -b buildtupledata.C+\(\"dtPbj60akPu4PF\"\) >out/dtPbj60akPu4PF &
root -l -q -b buildtupledata.C+\(\"dtPbj80akPu4PF\"\) >out/dtPbj80akPu4PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakPu4PF\"\) >out/dtPbbjtakPu4PF &
root -l -q -b buildtupledata.C+\(\"dtPbj60akPu3PF\"\) >out/dtPbj60akPu3PF &
root -l -q -b buildtupledata.C+\(\"dtPbj80akPu3PF\"\) >out/dtPbj80akPu3PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakPu3PF\"\) >out/dtPbbjtakPu3PF &

#root -l -q -b buildtupledata.C+\(\"dtPbj60akCs4PF\"\) >out/dtPbj60akCs4PF &
#root -l -q -b buildtupledata.C+\(\"dtPbbjtakCs4PF\"\) >out/dtPbbjtakCs4PF &
#root -l -q -b buildtupledata.C+\(\"dtPbj60akCs3PF\"\) >out/dtPbj60akCs3PF &
#root -l -q -b buildtupledata.C+\(\"dtPbbjtakCs3PF\"\) >out/dtPbbjtakCs3PF &

wait