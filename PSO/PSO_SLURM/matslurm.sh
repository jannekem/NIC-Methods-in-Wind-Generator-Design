#!/bin/sh
#SBATCH --time=0-00:05:00 --mem-per-cpu=500
#SBATCH -p play
#SBATCH -o job-%a.out
#SBATCH  --array=1-100

# define the working dir on the compute node
node_dir="local/$USER"

# check if it doesn't exist and create one
[ -d $node_dir ] || mkdir $node_dir

module load matlab
matlab -nojvm -r "PSO_function($SLURM_ARRAY_TASK_ID,50,500,50,0.7)"
