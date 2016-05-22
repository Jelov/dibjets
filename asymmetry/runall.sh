root -l <<< ".L tellmetruth.C+"

#root -l -b -q tellmetruth.C+\(\"results_0504_1_no_subtraction\",0,false,false\) > out_1 &
#root -l -b -q tellmetruth.C+\(\"results_0504_2_naive_subtraction\",1,false,false\) > out_2 &
#root -l -b -q tellmetruth.C+\(\"results_0504_3_estimated_subtraction\",2,false,false\) > out_3 &
#root -l -b -q tellmetruth.C+\(\"results_0504_4_trigger_corr\",2,true,false\) > out_4 &
#root -l -b -q tellmetruth.C+\(\"results_0504_5_tagging_corr\",2,true,true\) > out_5 &

#root -l -b -q tellmetruth.C+\(\"tmp\",2,false,false\)

#root -l -b -q tellmetruth.C+\(\"0520_dphi\"\)
root -l -b -q moneyplot.C+\(\"0520_dphi\"\)
