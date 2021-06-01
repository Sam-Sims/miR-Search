#!/bin/bash
#SBATCH --mail-user=naa17jnu@uea.ac.uk
#SBATCH --mail-type=END
cat accessions.txt | while read line
do
    echo "Reading ${line}"
    if [ ! -f /gpfs/home/naa17jnu/scratch/seq_runs/${line}.fastq.gz ] # check if file doesnt already exist
    then
        echo "File ${line} doesnt exist - downloading"
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/000/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 000"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/000/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/001/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 001"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/001/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 002"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/003/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 003"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/003/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/004/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 004"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/004/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/005/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 005"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/005/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/006/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 006"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/006/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/007/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 007"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/007/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
        if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/008/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 008"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/008/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
         if wget --spider ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/009/${line}/${line}.fastq.gz 2>/dev/null
        then
            echo "Found ${line}.fastq.gz on 008"
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/009/${line}/${line}.fastq.gz -P /gpfs/home/naa17jnu/scratch/seq_runs
        fi
    fi
done