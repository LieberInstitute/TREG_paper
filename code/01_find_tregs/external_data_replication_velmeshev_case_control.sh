#!/bin/bash
#$ -cwd
#$ -l mem_free=200G,h_vmem=200G,h_fsize=100G
#$ -N external_data_replication_velmeshev_case_control
#$ -o logs/external_data_replication_velmeshev_case_control.txt
#$ -e logs/external_data_replication_velmeshev_case_control.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript external_data_replication_velmeshev_case_control.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
