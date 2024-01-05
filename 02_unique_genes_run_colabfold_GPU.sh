#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job
#SBATCH --export=ALL
# Define job name
#SBATCH -J colab_GPU
# Define a standard output file. When the job is running, %N will be replaced by the name of
# the first node where the job runs, %j will be replaced by job id number.
#SBATCH -o colabfold_gpu.%N.%j.out
# Define a standard error file
#SBATCH -e colabfold_gpu.%N.%j.err
# Request the GPU partition (gpu). We don't recommend requesting multiple partitions, as the specifications of the nodes in these partitions are different.
#SBATCH -p gpu
# Request the number of nodes
#SBATCH -N 1
# Request the number of GPUs per node to be used (if more than 1 GPU per node is required, change 1 into Ngpu, where Ngpu=2,3,4)
#SBATCH --gres=gpu:1
# Request the number of CPU cores. (There are 24 CPU cores and 4 GPUs on each GPU node in partition gpu,
# so please request 6*Ngpu CPU cores, i.e., 6 CPU cores for 1 GPU, 12 CPU cores for 2 GPUs, and so on.)
#SBATCH -n 6
# Set time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores
##SBATCH --mem=90GB or #SBATCH --mem-per-cpu=9000M
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
##SBATCH --mail-user=username@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
#
# Set your maximum stack size to unlimited
ulimit -s unlimited
# Set OpenMP thread number
export OMP_NUM_THREADS=$SLURM_NTASKS

# Load tensorflow and relevant modules
module purge
module load apps/anaconda3/5.2.0
#use source activate gpu to get the gpu virtual environment
source activate gpu

# List all modules
module list

echo =========================================================
echo SLURM job: submitted  date = `date`
date_start=`date +%s`

hostname
echo Current directory: `pwd`

echo "Print the following environmetal variables:"
echo "Job name             : $SLURM_JOB_NAME"
echo "Job ID               : $SLURM_JOB_ID"
echo "Job user             : $SLURM_JOB_USER"
echo "CUDA_VISIBLE_DEVICES : $CUDA_VISIBLE_DEVICES"
echo "GPU_DEVICE_ORDINAL   : $GPU_DEVICE_ORDINAL"

echo "Running GPU jobs:"
echo Job output begins                                           
echo -----------------     

# This script was run on the Barkla cluster at Liverpool
# Files were stored on in /sharedscratch
# It should be run inside the fungiDB directory following on from the fasta processing step

# It runs colabfold for all fasta files in target directory
# 3 random seeds and includes relaxation step 
# For anyone wishing to run predictions on a local machine you can see script './code/run_colabfold_local.sh'
# If certain predictions fail and you wish to run colabfold for certain target predictions  you can see script './code/run_colabfold_cluster_GPU_FAILED_JOBS.sh'

#cd /mnt/lustre/users/ejohn10/fungiDB/unique_genes/
cd ./processed_fasta/unique_genes/

for i in *.fasta; do

  FASTA=$(basename ${i} .fasta)
  echo "Creating colabfold model for ${FASTA}.fasta"
  
  colabfold_batch --num-models 3 --use-gpu-relax --amber ${FASTA}.fasta ${FASTA}/ 
  # colabfold_batch --num-models 1 --use-gpu-relax --amber N1UWG6_L.intrerrogans_serovar_australis.fasta N1UWG6_L.intrerrogans_serovar_australis/ 

 done



#deactivate the gpu virtual environment
source deactivate gpu

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
