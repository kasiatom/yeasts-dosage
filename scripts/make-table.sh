#!/bin/bash

set -e -o pipefail
path="/home/kasia.tomala/clean-vep"

bcftools query -f \
'%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_Gene\t%INFO/VEP_SYMBOL\t%INFO/VEP_IMPACT\t%INFO/VEP_Consequence\t%INFO/VEP_BIOTYPE\t%INFO/VEP_Protein_position\t%INFO/VEP_Amino_acids\t%INFO/VEP_WORST_Gene\t%INFO/VEP_WORST_SYMBOL\t%INFO/VEP_WORST_IMPACT\t%INFO/VEP_WORST_Consequence[\t%GT]\n' \
"$path"/1011Matrix_vep_splited2.gvcf.gz \
-o "$path"/tmp.tsv 

printf "CHROM\tPOS\tREF\tALT\tVEP_Gene\tVEP_SYMBOL\tVEP_IMPACT\tVEP_Consequence\tVEP_BIOTYPE\tVEP_Protein_position\tVEP_Amino_acids\tVEP_WORST_Gene\tVEP_WORST_SYMBOL\tVEP_WORST_IMPACT\tVEP_WORST_Consequence\n" \
> "$path"/header1

bcftools view -h "$path"/1011Matrix_vep_splited2.gvcf.gz | tail -1 | cut -f 10- > "$path"/header2

cat <(paste "$path"/header1 "$path"/header2 ) "$path"/tmp.tsv > "$path"/1011Matrix_vep_splited2.tsv