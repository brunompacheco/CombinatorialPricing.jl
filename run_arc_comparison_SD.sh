#!/bin/bash
#SBATCH --mail-user=mpacheco.bruno@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-40
module load julia

FILES_EXPR="problems/interdiction-flatten/*.ki"

ls -1 $FILES_EXPR | shuf | awk -v count=$SLURM_ARRAY_TASK_COUNT -v id=$SLURM_ARRAY_TASK_ID 'NR % count == id - 1 {print}' | parallel julia --project arc_comparison.jl {}

echo "Job's done!"
