#!/bin/bash

## This script removes molecules with ambiguous BX codes, Phred mapping QUAL <50, and with less than four associated read pairs. 
## This script also adds sample ID (i.e., IID) to each molecule BX code for downstream analyses

##Example run
## ./filter_molecules.sh AQ15.P02.BX_sorted.linked_reads.full.bed

input=$1
output=`echo $input | cut -f1-3 -d"."`
IID=`echo $input | cut -f1 -d"."`

awk -F'\t' -v sample="${IID}_" 'BEGIN {OFS="\t"} {$4=sample $4; print}' $input > $output.labeled.bed
grep -vE "[ABD]00" $output.labeled.bed | awk '$5 >= 50 && $10 >= 4' | sort -k1,1 -k2,2n > $output.filter.bed
rm $output.labeled.bed
