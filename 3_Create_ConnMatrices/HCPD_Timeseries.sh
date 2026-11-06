#!/bin/bash

singularity run --cleanenv \
    insuladev/software/singularity_images/neuroimaging_packages_0.0.1.sif \
    Rscript --save insuladev/code/HCPD/3_Create_ConnMatrices/HCPD_Timeseries.R

    
