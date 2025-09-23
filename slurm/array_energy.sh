#!/usr/bin/env bash
#SBATCH -A PAS2277
#SBATCH --job-name=accpt_array
#SBATCH --output=out/corals_50km_10m_%a.out
#SBATCH --error=err/corals_50km_10m_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --array=1-7            # number of energy bins

cd ~/../../../fs/scratch/PAS2277/linton93/CoRaLS_MC/

#ALT=$1     # in km
#ICE=$2     # in m

# Slurm‐provided vars
IDX="$SLURM_ARRAY_TASK_ID"

# Pass task ID and task ID+1 to Julia
TASK1=$IDX

echo "[$SLURM_ARRAY_JOB_ID:$IDX]  alt=${ALT} km   ice=${ICE} m   task_range=${TASK1}-${TASK2}"
# call Julia with (altitude, ice_depth, bin_start, bin_end)
julia CoRaLS.jl/slurm/compute_acceptance.jl "$ALT" "$ICE" "$TASK1" 

#echo "Done!"
