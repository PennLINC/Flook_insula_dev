#!/bin/bash

#SBATCH --mem=64G\
#SBATCH    --output=/cbica/projects/insuladev/logs/subcortex/%x_%j.out \
#SBATCH    --error=/cbica/projects/insuladev/logs/subcortex/%x_%j.err \
#SBATCH    --job-name=Subcortex_Setup

insuladev/code/Subcortical/Subcortex_SetUp.sh