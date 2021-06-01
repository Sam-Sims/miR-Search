#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END
#SBATCH --mem=80G
#SBATCH -p compute-24-96

export PATH=/gpfs/home/naa17jnu/.conda/envs/r-env/bin:$PATH
source activate r-env
Rscript sleuth-new.R
