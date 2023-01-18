#!/usr/bin/bash
# WAJ 2021-01-12

# run on scg
#SBATCH --job-name=corrPlots
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=24:00:00
#SBATCH --mem=32G

# Loading Modules
module purge
ml R/4.1.2

# Log the output of the run
exec 1>> $(date -d "today" +"%Y%m%d%H%M")_command.log 2>&1
set -ex
start=`date +%s`

echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

# Run that bad boy!
Rscript --vanilla Correlation_Plots_and_Pathway_Enrichment.R

# Output Runtime
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
