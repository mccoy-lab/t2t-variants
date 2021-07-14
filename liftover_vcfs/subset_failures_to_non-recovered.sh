#!/bin/bash

recovered_file=$1
full_rejected_file=$2
not_recovered_output=$3

subset_script=subset_to_non-recovered.py

python python ${subset_script} ${recovered_file} ${full_rejected_file} ${not_recovered_output}
