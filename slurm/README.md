# slurm

This directory contains scripts for submitting jobs to a slurm batch system on an HPC (specifically the Ohio Supercomputer). The goal is to be able to easily submit jobs for running many simulations in parallel with long walltimes.

## Instructions

The main file is array_acceptance.sh . Submission to the slurm job manager looks like this: sbatch --array=<min>-<max> --export=ALL,ALT=<altitude>,ICE=<ice depth> array_acceptance.sh . 

The job script will run compute_acceptance.jl , which sets up a detector and reads in the variables from the job script to calculate the acceptance at a given energy bin. The energy bin is computed by multiplying 0.5 EeV by the number in the job array (<min> does not need to be 1, but it does need to be less than <max>). Currently, the width is linear, but in the future this may be updated to create half-decade energy bins.

### Example:

The below command will submit an array of 10 jobs, each one running simulations with a payload altitude of 50 km and an ice depth of 10 m. 
sbatch --array=1-10 --export=ALL,ALT=50,ICE=10 array_acceptance.sh

# TODO: 

Some work to consider adding:

1. Make energy bins half-decade in width
2. Accept ntrials as argument in job script
3. Generalize compute_acceptance.jl to accept more parameters (ex: orbital type, trigger, frequency range, etc). For now, those are easy to change in the script.
4. Make similar scripts for other quantities we may be interested in computing.
5. Add plotting scripts to read outputs from the array jobs and concatenate them into a single file and plot.
