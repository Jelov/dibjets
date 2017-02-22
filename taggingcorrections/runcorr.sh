root -l <<< ".L deriveOffTagEff.C+"

for CSV in 0.9 #0 0.1 0.2 #0.3 0.4 0.5 #0.6 0.7 0.8 #0 0.1 0.2 
do
  root -l -b -q deriveOffTagEff.C+\(false,false,true,false,0,$CSV\) &
  root -l -b -q deriveOffTagEff.C+\(true,false,true,false,0,$CSV\) &
  root -l -b -q deriveOffTagEff.C+\(true,false,true,false,1,$CSV\) &
  root -l -b -q deriveOffTagEff.C+\(true,false,true,false,2,$CSV\) &
  root -l -b -q deriveOffTagEff.C+\(true,false,true,false,3,$CSV\) &

  wait
done

