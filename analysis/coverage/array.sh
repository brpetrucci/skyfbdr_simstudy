#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=5-1:00:00 # probably overkill but 
#SBATCH --array=1-1000

#SBATCH --output=output/coverage/jobs/job_%A_%a.out
#SBATCH --error=output/coverage/jobs/job_%A_%a.err

#SBATCH --job-name="array_coverage"

#SBATCH --mail-user=petrucci@iastate.edu   # my e
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# move to the right directory
cd /work/LAS/phylo-lab/petrucci/skyfbdr_simstudy

# load modules
module load boost/1.81.0-m2umk6c

# source it, the parameter combination, and the actual script
rb analysis/coverage/master.Rev $SLURM_ARRAY_TASK_ID 
