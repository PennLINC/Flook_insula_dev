#!/bin/bash

singularity run --cleanenv \
    insuladev/software/singularity_images/modelarray_confixel_latest.sif \
    Rscript --save insuladev/code/HCPD/7_Change_1toNAN/change1toNAN.R