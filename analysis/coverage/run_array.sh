#!/bin/sh

for i in `seq 566 566`
do
printf "rep <- " >> analysis/coverage/aux_$i.Rev
echo $i >> analysis/coverage/aux_$i.Rev
rb analysis/coverage/aux_$i.Rev analysis/coverage/master.Rev
rm analysis/coverage/aux_$i.Rev
done
