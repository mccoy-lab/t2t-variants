#!/bin/sh

GATK=/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0/gatk
CHAIN=/scratch4/mschatz1/rmccoy22/CHM13/hg38.t2t-chm13-v1.0.over.chain
REF_FASTA=/scratch4/mschatz1/rmccoy22/CHM13/chm13.draft_v1.0.fasta

in_path=$1
lifted_out_path=$2
rejected_out_path=$3

$GATK LiftoverVcf -I ${in_path} -O ${lifted_out_path} --CHAIN ${CHAIN} -R ${REF_FASTA} --REJECT ${rejected_out_path}
