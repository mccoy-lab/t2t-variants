#!/bin/sh

GATK=/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0/gatk
CHAIN=/scratch4/mschatz1/rmccoy22/CHM13/hg38.t2t-chm13-v1.0.over.chain
REF_FASTA=/scratch4/mschatz1/rmccoy22/CHM13/chm13.draft_v1.0.fasta

in_path=$1
lifted_out_path=$2 # Should end with .gz
rejected_out_path=$3 #should end with .gz

$GATK LiftoverVcf -I ${in_path} -O ${lifted_out_path} --CHAIN ${CHAIN} -R ${REF_FASTA} --REJECT ${rejected_out_path}


# The below section may be necessary if LiftoverVcf encodes the output as ISO-8859-1,
# but should otherwise be unnecessary

gunzip -c ${rejected_out_path} > ${rejected_out_path}.nogz
iconv -f ISO-8859-1 -t UTF-8 ${rejected_out_path}.nogz -o ${rejected_out_path}.tmp
mv ${rejected_out_path}.tmp ${rejected_out_path}.nogz
gzip ${rejected_out_path}.nogz -c > ${rejected_out_path}
rm ${rejected_out_path}.nogz

gunzip -c ${lifted_out_path} > ${lifted_out_path}.nogz
iconv -f ISO-8859-1 -t UTF-8 ${lifted_out_path}.nogz -o ${lifted_out_path}.tmp
mv ${lifted_out_path}.tmp ${lifted_out_path}.nogz
gzip ${lifted_out_path}.nogz -c > ${lifted_out_path}
rm ${lifted_out_path}.nogz
