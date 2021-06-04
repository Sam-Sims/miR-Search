#!/bin/bash
cat input_mirs.txt | while read line
do
    echo "Reading ${line}"
    python mir-search ggplot-format -c ggplot_format_output_shape/${line}.csv RNA_fold_out/${line}.csv -o ${line}.csv
    echo "Done"
done