#!/bin/bash

singularity run --cleanenv \
    insuladev/software/singularity_images/modelarray_confixel_latest.sif \
    Rscript --save insuladev/code/HCPD/8_ROI_GAM/ROI_GAM_SetUp.R