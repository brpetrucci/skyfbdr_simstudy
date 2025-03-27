# go to the correct directory
cd .../skyfbdr_simstudy 

# get the rep value
rep=$SLURM_ARRAY_TASK_ID

# create a file to hold the rep
touch aux/aux_${1}_$rep.Rev

# echo the definitions on it
printf "rep <- " >> aux/aux_${1}_$rep.Rev
echo $rep >> aux/aux_${1}_$rep.Rev

# source it, the parameter combination, and the actual script
rb aux/aux_${1}_$rep.Rev simulation/accuracy/refs/ref_${1}.Rev analysis/accuracy/master.Rev

# remove file
rm aux/aux_${1}_$rep.Rev
