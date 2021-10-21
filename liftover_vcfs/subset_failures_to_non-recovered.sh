#!/bin/bash

recovered_file=$1
full_rejected_file=$2
not_recovered_output=$3 # should end with .gz

subset_script=subset_failures_to_non-recovered.py

python ${subset_script} ${recovered_file} ${full_rejected_file} ${not_recovered_output}.nogz

gzip -c ${not_recovered_output}.nogz > ${not_recovered_output}
rm ${not_recovered_output}.nogz
