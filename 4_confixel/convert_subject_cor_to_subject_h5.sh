#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=18G
#SBATCH --time=00:24:00
#SBATCH --output=insuladev/HCPDbold_1/subjectcsv_to_subjecth5_logs/convert-%A_%a.log
#SBATCH --array=1-563

micromamba activate modelarray
CSV=$(find insuladev/HCPDbold_1/cor_results -name '*.csv' | sort | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
python insuladev/code/HCPD/4_confixel/con_csv_to_h5.py $CSV


