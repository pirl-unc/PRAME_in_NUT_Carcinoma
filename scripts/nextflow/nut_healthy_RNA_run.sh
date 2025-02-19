#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=10g
#SBATCH --time=2-00:00 # time (D-HH:MM)
#SBATCH --job-name=KD_rna_all_reps_nextflow
#SBATCH --mail-user=email@unc.edu
#SBATCH --mail-type=end

## go to directory where you'd like temp files stored
cd /datastore/scratch/users/skpet/ 

# run the pipeline
nextflow run nf-core/rnaseq -r 3.14.0 -profile unc_lccc -params-file /scripts/nextflow/parameter_nut_healthy_lung_RNA.yaml

