#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python mir-search ggplot-format -s "sleuth_output_transcript/condition${line} .csv" icshape-align_out_flanking/${line}.pickle -o ${line}.csv
    echo "Done"
done