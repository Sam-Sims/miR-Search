#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python mir-search mir-search -m -i ${line} final_trans.fasta
    echo "Done"
done