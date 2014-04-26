#!/bin/sh
#SBATCH --time=0-02:00:00 --mem-per-cpu=500
#SBATCH -o job-%a.out
#SBATCH  --array=1-200

# change to working directory
cd $WRKDIR/PSO_SLURM

# launch computation
module load matlab
matlab -nojvm -r "PSO_function($SLURM_ARRAY_TASK_ID,50,500,50,0.7)"
