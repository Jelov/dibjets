#root -l <<< ".L buildtupledata.C+"
root -l <<< ".L buildtuplemc.C+"

#pp
root -l -q -b buildtuplemc.C+\(\"mcppqcdak4PF\"\) >out/mcppqcdak4PF &
root -l -q -b buildtuplemc.C+\(\"mcppbjtak4PF\"\) >out/mcppbjtak4PF &
root -l -q -b buildtuplemc.C+\(\"mcppqcdak3PF\"\) >out/mcppqcdak3PF &
root -l -q -b buildtuplemc.C+\(\"mcppbjtak3PF\"\) >out/mcppbjtak3PF &
#pyquen
#root -l -q -b buildtuplemc.C+\(\"mcPbpqcakPu4PF\"\) >out/mcPbpqcakPu4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbpfcakPu4PF\"\) >out/mcPbpfcakPu4PF &

#PbPb
root -l -q -b buildtuplemc.C+\(\"mcPbqcdakPu4PF\"\) >out/mcPbqcdakPu4PF &
root -l -q -b buildtuplemc.C+\(\"mcPbbjtakPu4PF\"\) >out/mcPbbjtakPu4PF &
root -l -q -b buildtuplemc.C+\(\"mcPbqcdakPu3PF\"\) >out/mcPbqcdakPu3PF &
root -l -q -b buildtuplemc.C+\(\"mcPbbjtakPu3PF\"\) >out/mcPbbjtakPu3PF &

#root -l -q -b buildtuplemc.C+\(\"mcPbqcdakPu3PF\"\) >out/mcPbqcdakPu3PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbjtakPu3PF\"\) >out/mcPbbjtakPu3PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbqcdakCs4PF\"\) >out/mcPbqcdakCs4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbjtakCs4PF\"\) >out/mcPbbjtakCs4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbqcdakCs3PF\"\) >out/mcPbqcdakCs3PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbjtakCs3PF\"\) >out/mcPbbjtakCs3PF &

# root -l -q -b buildtuplemc.C+\(\"mcPbqcdakPu4Calo\"\) >out/mcPbqcdakPu4Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbbjtakPu4Calo\"\) >out/mcPbbjtakPu4Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbqcdakPu3Calo\"\) >out/mcPbqcdakPu3Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbbjtakPu3Calo\"\) >out/mcPbbjtakPu3Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbqcdakCs4Calo\"\) >out/mcPbqcdakCs4Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbbjtakCs4Calo\"\) >out/mcPbbjtakCs4Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbqcdakCs3Calo\"\) >out/mcPbqcdakCs3Calo &
# root -l -q -b buildtuplemc.C+\(\"mcPbbjtakCs3Calo\"\) >out/mcPbbjtakCs3Calo &

 

wait