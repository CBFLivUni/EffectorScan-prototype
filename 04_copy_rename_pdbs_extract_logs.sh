#!/bin/bash

# Script to extract and rename PDB files 
# It also extracts the relevant log information and sequencing depth image
# It only extracts the PDB for the best ranked structure 
# Should be ran in the project directory

# Create directories
mkdir -p ./pdbs/unique_genes/sequencing_depth/
mkdir -p ./pdbs/unique_genes/logs/
mkdir -p ./pdbs/core_genes/sequencing_depth/
mkdir -p ./pdbs/core_genes/logs/

# Unique genes
for i in ./processed_fasta/unique_genes/*/ ; do
    
    NAME=$(basename ${i} /)
    echo "$NAME"
    cp ./processed_fasta/unique_genes/$NAME/*_relaxed_rank_001_*.pdb ./pdbs/unique_genes/$NAME.pdb
    cp ./processed_fasta/unique_genes/$NAME/*coverage.png ./pdbs/unique_genes/sequencing_depth/"$NAME"_coverage.png
    grep 'rank_001' ./processed_fasta/unique_genes/$NAME/log.txt  > ./pdbs/unique_genes/logs/$NAME.txt
    
done


# Core genes
for i in ./processed_fasta/core_genes/*/ ; do
    
    NAME=$(basename ${i} /)
    echo "$NAME"
    cp ./processed_fasta/core_genes/$NAME/*_relaxed_rank_001_*.pdb ./pdbs/core_genes/$NAME.pdb
    cp ./processed_fasta/core_genes/$NAME/*coverage.png ./pdbs/core_genes/sequencing_depth/"$NAME"_coverage.png
    grep 'rank_001' ./processed_fasta/core_genes/$NAME/log.txt  > ./pdbs/core_genes/logs/$NAME.txt
    
done

