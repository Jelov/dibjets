root -l <<< ".L tellmetruth.C+"


#root -l -b -q tellmetruth.C+\(\"0525_eclipseminus1GeV\"\)
#root -l -b -q tellmetruth.C+\(\"0527\"\)
root -l -b -q moneyplot.C+\(\"0527\"\)



#root -l -b -q tellmetruth.C+\(\"0524_BB\",true,true,true\) &
#root -l -b -q tellmetruth.C+\(\"0524_BX\",true,true,false\) &
#wait
#root -l -b -q tellmetruth.C+\(\"0524_XB\",true,false,true\) &
#root -l -b -q tellmetruth.C+\(\"0524_XX\",true,false,false\) &
#wait

# root -l -b -q tellmetruth.C+\(\"0524_BBnocorr\",false,true,true\) &
# root -l -b -q tellmetruth.C+\(\"0524_BXnocorr\",false,true,false\) &
# wait
# root -l -b -q tellmetruth.C+\(\"0524_XBnocorr\",false,false,true\) &
# root -l -b -q tellmetruth.C+\(\"0524_XXnocorr\",false,false,false\) &

