#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python mir-search icshape-align -s /home/sammy/Documents/MSc_Research_Project/9.transSHAPE/final.shape -t mir_search_output/${line}.csv -u final_trans.fasta --nopercent
    echo "Done"
done