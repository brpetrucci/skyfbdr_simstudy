#!/bin/sh

for j in `seq 63 63`
do
  for i in `seq 1 100`
  do
    start=`gdate +%s.%N`

    printf "rep <- " >> analysis/accuracy/aux_$i.Rev
    echo $i >> analysis/accuracy/aux_$i.Rev
    rb simulation/accuracy/refs/ref_$j.Rev analysis/accuracy/aux_$i.Rev analysis/accuracy/master.Rev
    rm analysis/accuracy/aux_$i.Rev
    
    end=`gdate +%s.%N`
    
    elapsed=$(echo "$end - $start" | bc)
    printf "%.9f\n" "$elapsed" >> simulation/accuracy/time_$j.tsv
  done
done
