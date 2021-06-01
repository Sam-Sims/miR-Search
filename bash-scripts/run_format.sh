#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python miR-Search.py format -i sleuth_output/sleuth_results_${line}.csv mir_search_output/${line}.csv -o ${line}.csv
    echo "Done"
done