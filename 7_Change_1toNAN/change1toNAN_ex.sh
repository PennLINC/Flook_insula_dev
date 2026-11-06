#!/bin/bash

#SBATCH --mem=16G\
#SBATCH    --output=insuladev/logs/NANs_Logs/%x_%j.out \
#SBATCH    --error=insuladev/logs/NANs_Logs/%x_%j.err \
#SBATCH    --job-name=NAN_job

insuladev/code/HCPD/7_Change_1toNAN/change1toNAN.sh