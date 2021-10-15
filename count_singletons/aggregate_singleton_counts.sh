#!/bin/bash

cd /scratch4/mschatz1/rmccoy22/rmccoy22/count_singletons

for i in {1..22}; do cat chr${i}_CHM13_singletons.txt ; done >> CHM13_singletons.txt

head -1 chr1.singletons > 1KGP_singletons.txt
for i in {1..22}; do cat chr${i}.singletons | sed -e1,1d ; done >> 1KGP_singletons.txt

wc -l CHM13_singletons.txt | awk '{print $1"\tCHM13"}' > singleton_total.txt
sed -e1,1d 1KGP_singletons.txt | awk '{if ($3 == "S") print $5}' | sort | uniq -c | awk '{print $1"\t"$2}' >> singleton_total.txt

