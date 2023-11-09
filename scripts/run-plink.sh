#!/bin/bash

set -e -o pipefail

cd $HOME/distances

##define sample names
samples=$(zcat $HOME/clean-vep/1011Matrix_vep_splited2.gvcf.gz | grep -m1 "^#C" | sed 's/\t/ /g' )

## change X to XVII to cheat plink (X is not sex chromosome in yeasts)
zcat $HOME/clean-vep/1011Matrix_vep_splited2.gvcf.gz | sed 's/^X\t/XVII\t/' | bgzip > tmp.vcf.gz

## run plink2 --sample-diff
./plink2 \
        --vcf tmp.vcf.gz \
        --out 1011 \
        --allow-extra-chr \
        --sample-diff  \
         ids=$samples

## Hamming distance plink1.9
plink \
        --vcf tmp.vcf.gz \
        --out 1011-dist \
        --allow-extra-chr \
        --distance 1-ibs allele-ct square

## reformat mdist output (1-ibs) => from matrics to table
## ids
for i in {1..301};do cat 1011-dist.mdist.id >> ids2.txt;done
while read f; do  for i in {1..301}; do echo $f >> ids1.txt;done;done<1011-dist.mdist.id

## results
for i in {1..301}; do cut -f $i 1011-dist.mdist >> mist.data;done 

## merge
paste <(cut -f1 -d ' ' ids1.txt) <(cut -f1 ids2.txt) mist.data > 1011-ibs-distance.txt

## King table - with IBS0, HETHET and KINSHIP metrics 
./plink2 \
        --vcf tmp.vcf.gz \
        --allow-extra-chr \
        --out 1011-king \
        --make-king-table

## clean
rm tmp.cf.gz ids1.txt ids2.txt mist.data

