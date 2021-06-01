#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END

#Before running generate index file - kallisto index -i filename

# DO module load kallisto/0.46.1
module load kallisto/0.46.1
line=$1
echo "Reading ${line}"
if [ -f /gpfs/home/naa17jnu/scratch/fastp/${line}_clean.fastq.gz ] # check clean file actually exists
then
echo "File exists at /gpfs/home/naa17jnu/scratch/fastp/${line}_clean.fastq.gz" 
    if [ ! -f ${line}_block ] # check kallisto isnt in progress
    then
        echo "Running kallisto on ${line}"
        touch ${line}_block
        kallisto quant -i /gpfs/home/naa17jnu/scratch/kallisto/homo_sapiens_GRCh38 --single --fragment-length=180 --sd=20 --bias -b 100 -o /gpfs/home/naa17jnu/scratch/kallisto/${line}_kallisto /gpfs/home/naa17jnu/scratch/fastp/${line}_clean.fastq.gz 2>&1 | tee ${line}_stdout.txt
    else
        echo "Block file found - is kallisto running already?"
    fi
else
    echo "Cleaned file not found"
fi
rm -f ${line}_block
done
