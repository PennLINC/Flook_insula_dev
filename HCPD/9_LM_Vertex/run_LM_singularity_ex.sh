#!/bin/bash

#SBATCH --job-name=LM_Array_Job
#SBATCH --output=insuladev/logs/LM_Logs/%A_%a.out 
#SBATCH --error=insuladev/logs/LM_Logs/%A_%a.err
#SBATCH --mem=16G
#SBATCH --array=1-1478

insuladev/code/HCPD/9_LM_Vertex/run_LM_singularity.sh $SLURM_ARRAY_TASK_ID
