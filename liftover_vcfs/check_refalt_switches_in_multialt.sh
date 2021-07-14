#!/bin/bash

parser="parse_mult_alt_alleles.py"
resolver="mult_alt_resolver.py"
GATK="/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0/gatk"
CHAIN="/scratch4/mschatz1/rmccoy22/CHM13/hg38.t2t-chm13-v1.0.over.chain"
REF="/scratch4/mschatz1/rmccoy22/CHM13/chm13.draft_v1.0.fasta"

foi=$1
parsed_multiallelic=$1.parsed
tmp_out_lifted=$1.parsed.lifted
tmp_out_rejected=$1.parsed.rejected
recovered_info=$1.parsed.recovered.info
recovered_out_file=$2

zcat ${foi} | python ${parser} | gzip -c > ${parsed_multiallelic}
${GATK} LiftoverVcf -I ${parsed_multiallelic} -O ${tmp_out_lifted} --CHAIN ${CHAIN} -R ${REF} --REJECT ${tmp_out_rejected} --RECOVER_SWAPPED_REF_ALT --WRITE_ORIGINAL_POSITION --WRITE_ORIGINAL_ALLELES

zcat ${tmp_out_lifted} | grep ^[^#] | cut -f 1-5,8 > ${recovered_info}

python ${resolver} ${recovered_info} ${foi} ${tmp_out_lifted} ${recovered_out_file}

rm ${parsed_multiallelic}
rm ${tmp_out_lifted}
rm ${tmp_out_rejected}
rm "${tmp_out_lifted}.tbi"
rm ${recovered_info}