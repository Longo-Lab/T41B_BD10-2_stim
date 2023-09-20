#!/usr/bin/env bash
# WAJ 2021-01-21

#SBATCH --job-name=DTUrtle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=batch
#SBATCH --account=flongo
#SBATCH --time=24:00:00
#SBATCH --mem=128G

module purge
ml R/4.1.2

# log the output of the run
exec 1>> $(date -d "today" +"%Y%m%d%H%M")_command.log 2>&1
set -ex

Rscript --vanilla DTUrtle_pipeline.R \
  --run_dge \
  --check_bias \
  -c genotype \
  -g WT \
  -s noSTIM \
  -j 32

squeue -j $SLURM_JOB_ID --Format=TimeUsed
