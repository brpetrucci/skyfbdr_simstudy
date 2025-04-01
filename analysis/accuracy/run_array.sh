#!/bin/sh

for i in `seq 1 1`
do
printf "rep <- " >> analysis/accuracy/aux_$i.Rev
echo $i >> analysis/accuracy/aux_$i.Rev
rb simulation/accuracy/refs/ref_62.Rev analysis/accuracy/aux_$i.Rev analysis/accuracy/master.Rev
rm analysis/accuracy/aux_$i.Rev
done
