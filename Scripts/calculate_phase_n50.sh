#!/bin/bash
#usage - ./calculate_phase_n50.sh <input example.n50.bed>
#conda activate base_genomics

file=$1
runID=`basename $file .bed`
samid=`echo $runID | sed 's/_.*//' | cut -f1 -d'.'`

# Read data from column 6 of example.n50.bed and print the lowest value

n50=($(sort -k 6,6n $file | awk '{print $6}' | head -n 1))

# Print the N50 value
echo "$samid N50 haplotype length is $n50 bases"