#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END

cat accessions.txt | while read line
do
    echo "Cleaning ${line}"
    ./fastp -i /gpfs/home/naa17jnu/scratch/seq_runs/${line}.fastq.gz -o /gpfs/home/naa17jnu/scratch/fastp/${line}_clean.fastq.gz -Q -L -j ${line}.json - h ${line}.html -R "${line}"
    echo "Done"
done