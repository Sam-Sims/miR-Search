#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END


# MUST DO module load fastqc/0.11.9 before running
module load fastqc/0.11.9
cat accessions.txt | while read line
do
    echo "Cleaning ${line}"
    fastqc /gpfs/home/naa17jnu/scratch/seq_runs/${line}.fastq.gz --outdir /gpfs/home/naa17jnu/scratch/fastqc_analysis
    echo "Done"
done