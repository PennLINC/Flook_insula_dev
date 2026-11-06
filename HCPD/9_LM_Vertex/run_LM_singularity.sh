#!/bin/bash

voxel_id=$1
singularity run --cleanenv \
    insuladev/software/singularity_images/modelarray_confixel_latest.sif \
    Rscript --save insuladev/code/HCPD/9_LM_Vertex/LM_array.R $voxel_id
    
