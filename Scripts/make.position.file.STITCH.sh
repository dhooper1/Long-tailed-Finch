#!/bin/bash

chrom=$1

#Extract SNP information for variants above SNP QUAL threshold >=500

bcftools query -f '%CHROM %POS %REF %ALT %QUAL\n' $chrom.concatenated.biallelic.repeatmask.vcf.gz | awk '$5>=500' | awk '{print $1,$2,$3,$4}' > $chrom.Q500.repeatmask.pos.txt;

sed -i -- $'s/ /\t/g' $chrom.Q500.repeatmask.pos.txt;

wc -l $chrom.Q500.repeatmask.pos.txt;

echo "Finished extracting high-quality SNP information from $chrom";