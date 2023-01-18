#!/usr/bin/bash
# WAJ 2021-01-14

# run on sherlock
#SBATCH --job-name=DTU_pre1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=24:00:00
#SBATCH --mem=128G

module purge
ml anaconda && source activate DTUrtle

parallel -j 16 cp /oak/stanford/scg/lab_flongo/t41b_BD10-2_stim/T41B_BD10_2_RNAseq/BD10_2_{}/*trimmed.fastq.gz /oak/stanford/scg/lab_flongo/t41b_BD10-2_stim/DTUrtle/fastq/ ::: '1' '2' '3' '4' '9' '10' '11' '12' '17' '18' '19' '20' '25' '26' '27' '28'

squeue -j $SLURM_JOBID --Format=TimeUsed
