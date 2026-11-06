#!/bin/bash

#SBATCH --mem=32G\
#SBATCH    --output=insuladev/logs/ROI_GAM/%x_%j.out \
#SBATCH    --error=insuladev/logs/ROI_GAM/%x_%j.err \
#SBATCH    --job-name=ROI_setup_job

insuladev/code/HCPD/8_ROI_GAM/ROI_GAM_SetUp.sh