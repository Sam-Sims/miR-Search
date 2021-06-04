#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python mir-search rnafold -t mir_search_output/${line}.csv -u final_trans.fasta -s "sleuth_output_transcript/condition${line} .csv" -f 30 -o ${line}.csv
    echo "Done"
done