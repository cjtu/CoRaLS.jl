#!/usr/bin/env bash
#SBATCH -A PAS0654
#SBATCH --job-name=accpt_array
#SBATCH --output=out/corals_%A_%akm_6m.out
#SBATCH --error=err/corals_%A_%akm_6m.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=4G

## Example: sbatch --array=1-20 --export=ALL,ENERGY=1,ICE=5 array_altitude.sh
## Example: sbatch --array=1-20 --export=ALL,ALT=1,ENERGY=1,ICE=5,ANT=4,TRIG=4,ANG=-90,FREQ1=300,TEXP=7,VAR=ALT array_altitude.sh

cd /users/PAS0654/machtay1/new_corals/

## Set the variables we read in

#ALT=$1					# altitude in km
#ENERGY=$2			# Energy multiplier (set to 1)
#ICE=$3					# ice depth in m
#ANT=$4					# number of antennas
#TRIG=$5				# number of triggers required
#ANG=$6					# pointing angle of antennas in degrees
#FREQ1=$7				# minimum frequency in band in MHz
#TEXP=$8				# exponent on 10 for ntrials


# Slurm‐provided vars
TASK="$SLURM_ARRAY_TASK_ID"

## Whichever of the above variables you want to vary with task
##  set equal to $TASK here:
eval ${VAR}=${TASK}
echo $VAR
echo ${ALT}

echo "alt=$((5 * ${ALT})) km   energyMult=${ENERGY}   ice=${ICE} m   ant=${ANT}   trig=${TRIG}   angle=${ANG} deg  freqMin=${FREQ1} MHz   TEXP=${TEXP}"
# call Julia with (altitude, ice_depth, bin_start, bin_end)
julia CoRaLS.jl/slurm/generic_acceptance.jl "$ALT" "$ENERGY" "$ICE" "$ANT" "$TRIG" "$ANG" "$FREQ1" "$TEXP" 

#echo "Done!"
