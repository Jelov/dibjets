#Build FCR (and BFA)
root -l <<< ".L buildtuplemc.C+"

root -l -q -b buildtuplemc.C+\(\"mcp1bfcak4PF\"\) >out/mcp1bfcak4PF &
root -l -q -b buildtuplemc.C+\(\"mcp2bfcak4PF\"\) >out/mcp2bfcak4PF &
root -l -q -b buildtuplemc.C+\(\"mcp3bfcak4PF\"\) >out/mcp3bfcak4PF &

#root -l -q -b buildtuplemc.C+\(\"mcppbfcak4PF\"\) >out/mcppbfcak4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu4PF\"\) >out/mcPbbfcakPu4PF &
# root -l -q -b buildtuplemc.C+\(\"mcppbfcak4PF\"\) >out/mcppbfcak4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu4PF\"\) >out/mcPbbfcakPu4PF &

wait

#root -l -q -b buildtuplemc.C+\(\"mXppbfcak4PF\"\) >out/mcppbfcak4PF &
#root -l -q -b buildtuplemc.C+\(\"mXPbbfcakPu4PF\"\) >out/mXPbbfcakPu4PF &


#root -l -q -b buildtuplemc.C+\(\"mcppbfcak3PF\"\) >out/mcppbfcak3PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu3PF\"\) >out/mcPbbfcakPu3PF &



#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs4PF\"\) >out/mcPbbfcakCs4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs3PF\"\) >out/mcPbbfcakCs3PF &

#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu4Calo\"\) >out/mcPbbfcakPu4Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu3Calo\"\) >out/mcPbbfcakPu3Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs4Calo\"\) >out/mcPbbfcakCs4Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs3Calo\"\) >out/mcPbbfcakCs3Calo &


wait