#!/bin/bash

export TF_FORCE_UNIFIED_MEMORY="1"
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
export TF_FORCE_GPU_ALLOW_GROWTH="true"

# Run colabfold for all fasta files in directory
# 3 random seeds and includes relaxation step 
for FASTA in $(cat unique_genes_failed_jobs.txt); do


  echo "Creating colabfold model for ${FASTA}.fasta"
  
  colabfold_batch --num-models 3 --use-gpu-relax --amber ${FASTA}.fasta ${FASTA}/ 
  # colabfold_batch --num-models 1 --use-gpu-relax --amber N1UWG6_L.intrerrogans_serovar_australis.fasta N1UWG6_L.intrerrogans_serovar_australis/ 
  

done