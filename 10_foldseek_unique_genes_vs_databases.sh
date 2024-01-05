#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J foldseek
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o foldseek.%u.%N.%j.out
# Define a standard error file
#SBATCH -e foldseek.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 20
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores 
##SBATCH --mem-per-cpu=9000M
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
##SBATCH --mail-user=username@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
#
# Set the maximum stack size to unlimited
ulimit -s unlimited
# Set OpenMP thread number
export OMP_NUM_THREADS=$SLURM_NTASKS

# Load your own modules
module purge

# List all modules
module list
#
#
echo =========================================================   
echo SLURM job: submitted  date = `date`      
date_start=`date +%s`

echo -------------  
echo Job output begins                                           
echo -----------------                                           
echo

hostname

echo "Print the following environmetal variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job user                     : $SLURM_JOB_USER"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"
echo "Submit directory             : $SLURM_SUBMIT_DIR"
echo "Temporary directory          : $TMPDIR"
echo "Submit host                  : $SLURM_SUBMIT_HOST"
echo "Queue/Partition name         : $SLURM_JOB_PARTITION"
echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Hostname of 1st node         : $HOSTNAME"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

echo   
echo "Running job:"
echo   

echo =========================================================   


# Script to query unique effector proteins against various protein DBs using FoldSeek 
#foldseek easy-search --taxon-list "4751" FOXG_02842.pdb databases/swissprot fungi_test tmp --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" 


# Create directories
mkdir -p outputs/structural_alignment/all_unique_genes_vs_swissprot
mkdir -p outputs/structural_alignment/all_unique_genes_vs_pdb
mkdir -p outputs/structural_alignment/all_unique_genes_vs_proteome
mkdir -p outputs/structural_alignment/all_unique_genes_vs_uniprot50

# Swiss-prot
for ID in $(cat outputs/unique_genes_IDs.txt); do

  echo "Querying ${ID}.pdb against SwissProt database for structural homology"
  foldseek easy-search \
  --threads $SLURM_NTASKS \
  --taxon-list "4751" \
  ./pdbs/unique_genes/$ID.pdb \
  databases/swissprot \
  swissprot_structural_homology_$ID.m8 tmp \
  --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" \
  --exhaustive-search

  mv *structural_homology*.m8 ./outputs/structural_alignment/all_unique_genes_vs_swissprot

done


 # PDB
for ID in $(cat outputs/unique_genes_IDs.txt); do

  echo "Querying ${ID}.pdb against PDB database for structural homology"
  foldseek easy-search \
  --threads $SLURM_NTASKS \
  --taxon-list "4751" \
  ./pdbs/unique_genes/$ID.pdb \
  databases/pdb \
  pdb_structural_homology_$ID.m8 tmp \
  --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" \
  --exhaustive-search

  mv *structural_homology*.m8 ./outputs/structural_alignment/all_unique_genes_vs_pdb

done


# Proteome
for ID in $(cat outputs/unique_genes_IDs.txt); do

  echo "Querying ${ID}.pdb against Proteome database for structural homology"
  foldseek easy-search \
  --threads $SLURM_NTASKS \
  --taxon-list "4751" \
  ./pdbs/unique_genes/$ID.pdb \
  databases/proteome \
  proteome_structural_homology_$ID.m8 tmp \
  --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" \
  --exhaustive-search

  mv *structural_homology*.m8 ./outputs/structural_alignment/all_unique_genes_vs_proteome

done

# Uniprot50
for ID in $(cat outputs/unique_genes_IDs.txt); do

  echo "Querying ${ID}.pdb against Uniprot50 database for structural homology"
  foldseek easy-search \
  --threads $SLURM_NTASKS \
  --taxon-list "4751" \
  ./pdbs/unique_genes/$ID.pdb \
  databases/uniprot50 \
  uniprot50_structural_homology_$ID.m8 tmp \
  --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" #\
  #--exhaustive-search , exhaustive search isn't used here due to size of DB

  mv *structural_homology*.m8 ./outputs/structural_alignment/all_unique_genes_vs_uniprot50

done


echo =========================================================   
# the ret flag is the return code, so you can spot easily if your code failed.
ret=$?

echo   
echo ---------------                                           
echo Job output ends                                           
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================   
echo SLURM job: finished   date = `date`   
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================   
exit $ret