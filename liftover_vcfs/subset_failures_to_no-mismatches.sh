#!/bin/bash

rejected_file=$1
output_position_failed=$2 # should end with .gz

subset_script=subset_failures_to_no_mismatches.py

python ${subset_script} ${rejected_file} | gzip -c > ${output_position_failed}
