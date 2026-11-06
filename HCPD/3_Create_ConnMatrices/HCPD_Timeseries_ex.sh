#!/bin/bash

qsub -l h_vmem=16G,s_vmem=16G \
    -b y \
    -V \
    -j n \
    -o insuladev/logs/ConMatrix_Logs \
    -e insuladev/logs/ConMatrix_Logs \
   insuladev/code/HCPD/3_Create_ConnMatrices/HCPD_Timeseries.sh

    