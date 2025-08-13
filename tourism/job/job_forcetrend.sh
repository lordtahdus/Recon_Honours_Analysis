#!/bin/bash
#SBATCH --job-name=sim_chunk
#SBATCH --output=tourism/job/logs/output_%A_%a.log
#SBATCH --error=tourism/job/logs/error_%A_%a.log
#SBATCH --array=98-133
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# exclude 19th iter since error in last run results
# Load the R module
module load R/4.4.0-mkl

export R_LIBS_USER=~/yi61/tsuu0007/R/library

# Calculate simulation index
INDEX=$(( $SLURM_ARRAY_TASK_ID ))

# Run the script
# original run
# Rscript ~/yi61/tsuu0007/Recon_Honours_Analysis/tourism/job/tourism_run.R $INDEX

# forcetrend run
Rscript ~/yi61/tsuu0007/Recon_Honours_Analysis/tourism/job/tourism_run_forcetrend.R $INDEX