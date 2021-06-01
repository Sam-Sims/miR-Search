# miR-Search

mir-search is a pipeline devloped to examine and analyse the role of miRNA target site structure in the context of gene expression

# Introduction

mir-search is a pipeline devloped as part of my MSc research project. It combines RNA-Seq data and [SHAPE-Seq](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/shape-seq.html) data in order to visualise the role that mRNA secondary structure has on miRNA mediated gene regulation.

mir-search is written in Python 3.8.

The mir-search pipeline consists of 4 modules:
- pymart
- mir-search
- icshape-align
- ggplot-format

#### pymart:
The pymart module allows you to download 3'UTR sequence data directly from biomart for a given subset of genes. It also cleans the 3'UTR data to remove unavailible sequence data and prepares it for direct use within mir-search module.

#### mir-search:
The mir-search module is used to identify all 6mer, 7mer and 8mer targets of a given micro-RNA

#### icshape-align:
The icshape-align module takes the output from the mir-search module as well as SHAPE-Seq data and extracts structure scores for each nucleotide in the target region

#### ggplot-format:
ggplot-format is a helper module that takes the output from the mir-search or icshape-align module and parses the data in a way that allows for direct plotting with the R library, ggplot2

#### rna-fold:
rna-fold is a module that implements the ViennaRNA Package to calculate the minimum free energy (MFE) structure for target site specified in the mir-search output. This takes the target site locations +/- a flanking region and pipes them to the RNAfold stdin, runs RNAfold and then captures the stdout. 

Also prints the MFE structure in dot-bracket notation and its free energy to stdout as well as the partition function and the base pairing probability matrix.
# Usage

#### pymart:
    
Usage: `python mir-search pymart -i "genelist" - o "output" -a`

Note: The genelist should be a csv file with the first column the ensembl gene IDs

template.xml contains the xml request string. This currently obtains 3'UTR sequences (specified via: `<Attribute name = "3utr" />`).

#### mir-search:

Usage: `python mir-search mir-search -i -m input_mir_name input_utr.fasta`

Note: The -m flag specifies the input micro RNA name as seen on mirbank.

#### icshape-align

Usage: `mir-search icshape-align -s icshape_output -t mir_search_output -u UTR file -f 0`

Note: The -f flag specifies the length of the flanking region to read