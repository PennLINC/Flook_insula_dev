#!/bin/bash
#manually update to run Subcortex_Setup.R, Subcortex_Setup2.R, Subcortex_Setup3.R, Subcortex_Setup4.R, Subcortex_Setup5.R

singularity run --cleanenv \
    insuladev/software/singularity_images/modelarray_confixel_latest.sif \
    Rscript --save insuladev/code/Subcortical/Subcortex_Setup3.R