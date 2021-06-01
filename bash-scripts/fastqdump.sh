#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END

source activate icshape-pipe
fasterq-dump SRR11164866 -O /gpfs/home/naa17jnu/scratch/fastq_output
