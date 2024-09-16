#!/bin/bash

## Calculate molecule N50 from output of bed_write.pl script after filtering results with filter_molecules.sh

input=$1
sample=`echo $input | cut -f1 -d"_"`

awk '{print $3 - $2}' $input | sort -nr > $sample.mol_lengths.tmp

total_sum=$(awk '{sum += $1} END {print sum}' $sample.mol_lengths.tmp)
half_sum=$(echo "$total_sum / 2" | bc -l)

n50_value=$(awk -v half_sum="$half_sum" '{cumsum += $1; if (cumsum >= half_sum) {print $1; exit}}' $sample.mol_lengths.tmp)
awk -v half_sum="$half_sum" '{cumsum += $1; if (cumsum >= half_sum) {print $1; exit}}' $sample.mol_lengths.tmp > $sample.mol.N50.out

rm $sample.mol_lengths.tmp

echo "$sample molecule N50 is: $n50_value"