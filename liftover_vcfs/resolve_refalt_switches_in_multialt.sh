#!/bin/bash

parser="parse_mult_alt_alleles.py"
resolver="mult_alt_resolver.py"
GATK="/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0/gatk"
CHAIN="/scratch4/mschatz1/rmccoy22/CHM13/hg38.t2t-chm13-v1.0.over.chain"
REF="/scratch4/mschatz1/rmccoy22/CHM13/chm13.draft_v1.0.fasta"

foi=$1 # This should be the rejected vcf file from the initial liftover
recovered_out_file=$2 # should end with .gz

parsed_multiallelic=$1.parsed.vcf.gz
tmp_out_lifted=$1.parsed.lifted.vcf.gz
tmp_out_rejected=$1.parsed.rejected.vcf.gz
recovered_info=$1.parsed.recovered.info


zcat ${foi} | python ${parser} | gzip -c > ${parsed_multiallelic}
${GATK} LiftoverVcf -I ${parsed_multiallelic} -O ${tmp_out_lifted} --CHAIN ${CHAIN} -R ${REF} --REJECT ${tmp_out_rejected} --RECOVER_SWAPPED_REF_ALT --WRITE_ORIGINAL_POSITION --WRITE_ORIGINAL_ALLELES

zcat ${tmp_out_lifted} | iconv -f iso8859-1 -t utf8 | gzip -c > ${tmp_out_lifted}.corrected
mv ${tmp_out_lifted}.corrected ${tmp_out_lifted}

zcat ${tmp_out_lifted} | grep -a ^[^#] | cut -f 1-5,8 > ${recovered_info}

python ${resolver} ${recovered_info} ${foi} ${tmp_out_lifted} ${recovered_out_file}.nogz

gzip ${recovered_out_file}.nogz -c > ${recovered_out_file}

rm ${parsed_multiallelic}
rm ${tmp_out_lifted}
rm ${tmp_out_rejected}
rm "${tmp_out_lifted}.tbi"
rm ${recovered_info}
rm ${recovered_out_file}.nogz
