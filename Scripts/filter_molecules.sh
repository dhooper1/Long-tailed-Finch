#!/bin/bash

## This script removes molecules with ambiguous BX codes, Phred mapping QUAL <50, and with less than four associated read pairs

input=$1
output=`echo $input | cut -f1-3 -d"."`

awk '!/A00/ && !/B00/ && !/D00/ && $5 >= 50 && $10 >= 4' $input > $output.filter.bed
