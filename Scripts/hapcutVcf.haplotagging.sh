#!/bin/bash

#This script is used to create a combined, phased VCF file for a single individual sample
#Script written by YFC and modified by DMH
#It does so with the following steps:
#
#1. generate two temporary VCF files - one with only heterozygous sites, and the other one with all sites
#2. HAPCUT2 pipeline, with basically three steps:
#2a. Parse the bam file for reads that correspond to the heterozygous sites for HAPCUT2 using the extractHAIRs utility with the 10X flag on to generate a BX-tagged, unlinked fragment file
#2b. Link the fragments using HAPCUT2's LinkFragments.py utility.
#2c. Run HAPCUT2 proper, with VCF output option
#3. Perform some basic data extraction from the HAPCUT2 output
#4. Merge the full VCF file with the HAPCUT2-derived VCF files and annotating the INFO and FORMAT fields accordingly
#5. Clean-up

bam=$1
vcf=$2
region=$3
runBC=`basename $bam .bam`
samid=`echo $runBC | sed 's/_.*//' | cut -f1 -d'.'`
sampleHetVCF=$samid.pMarkdup.$region.PL.AD.het.vcf
sampleVCF=$samid.pMarkdup.$region.PL.AD.vcf
chrom=`echo $region | sed 's/:/\t/' | cut -f 1`;
bin_dir=/mendel-nas1/dhooper/bin/miniconda3/envs/base_genomics/bin
out_dir=/mendel-nas1/dhooper/SNPs/phased_vcfs/$chrom
tmp=/mendel-nas1/dhooper/SNPs/phased_vcfs/fragments/$samid
if [ ! -e "$out_dir" ]; then mkdir -p $out_dir; fi
if [ ! -e "$tmp" ]; then mkdir -p $tmp; fi

echo "BAM: $bam"
echo "vcf: $vcf"
echo "runBC: $runBC"
echo "samID: $samid"
echo "chrom: $chrom"
echo "region: $region"
echo "sampleHetVCF: $sampleHetVCF"
echo "sampleVCF: $sampleVCF"
rm $tmp/$samid.$region.*

#1. generate two temporary VCF files - one with only heterozygous sites, and the other one with all sites
echo "Step 1: Writing temporary VCF files:"
echo "Step 1a: Write sample VCF with only heterozygous sites..."
$bin_dir/bcftools view -s $samid -r $region $vcf | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > $tmp/$sampleHetVCF

echo "Step 1b: Write sample VCF with all sites..."
$bin_dir/bcftools view -s $samid -r $region $vcf | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0}; !/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; !/^#/ && $10~/^0\/1/ {$9=substr($9, 4); $10=substr($10,5);print $0}' > $tmp/$sampleVCF

#Option to use if PG field has not been used as backup GT
##/fml/chones/local/bin/bcftools view -s $samid $vcf -i 'INFO/INFO_SCORE >= 0.2' | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[2] > GP[1] && GP[2] > GP[3]) {$10="0/1"substr($10,4);print $0}}' > $tmp/$sampleHetVCF
##/fml/chones/local/bin/bcftools view -s $samid $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10=$10":0|0";$9=$9":PG";print $0};  !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {$10=$10":1|1";$9=$9":PG";print $0}; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[1] > GP[2] && GP[1] > GP[3]) {$10=$10":0|0";$9=$9":PG";print $0}; if(GP[3]>GP[2] && GP[3] > GP[1]) {$10=$10":1|1";$9=$9":PG";print $0};if (GP[2] > GP[1] && GP[2] > GP[3]) {print $0}; if (GP[1] == GP[2] || GP[1] == GP[3] || GP[2] == GP[3]){print $0} }' > $tmp/$sampleVCF

#2. HAPCUT2 pipeline, with basically three steps:
echo "Step2: Running HAPCUT2:"

#2a. Parse the bam file for reads that correspond to the heterozygous sites for HAPCUT2 using the extractHAIRs utility with the 10X flag on to generate a BX-tagged, unlinked fragment file
echo "Step 2a: Extracting reads from BAM file that overlap het sites..."
##Currently uses conda build of hapcut2
extractHAIRS --10X 1 --bam $bam --VCF $tmp/$sampleHetVCF --region $region --out $tmp/$samid.$region.unlinked.fragments
grep -ax '.*' $tmp/$samid.$region.unlinked.fragments > $tmp/$samid.$region.unlinked.fragments.1
mv $tmp/$samid.$region.unlinked.fragments.1 $tmp/$samid.$region.unlinked.fragments

#2b. Link the fragments using HAPCUT2's LinkFragments.py utility
echo "Step 2b: Linking fragments using HAPCUT2's LinkFragments.py utility..."
##Currently uses conda python env and local build of hapcut2
pythonPATH=/mendel-nas1/dhooper/bin/miniconda3/envs/python/bin/
$pythonPATH/python /mendel-nas1/dhooper/bin/HapCUT2/utilities/LinkFragments.py --bam $bam --VCF $tmp/$sampleHetVCF --fragments $tmp/$samid.$region.unlinked.fragments --out $tmp/$samid.$region.linked.fragments -d 50000
#python3 LinkFragments.py --bam $bam --VCF $tmp/$sampleHetVCF --fragments $tmp/$out.$region.unlinked.fragments --out $tmp/$out.$region.linked.fragments -d 50000;

#2c. Run HAPCUT2 proper, with VCF output option
echo "Step 2c: Running HAPCUT2 proper, with VCF output..."
##Currently uses conda build of hapcut2
HAPCUT2 --fragments $tmp/$samid.$region.linked.fragments --VCF $tmp/$sampleHetVCF --out $tmp/$samid.$region.hapcut2.threshold_30.output --nf 1 --threshold 30 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1

##
echo "Processing HAPCUT2 output..."
/mendel-nas1/dhooper/SNPs/phased_vcfs/hapcut2_to_bed.sh $tmp/$samid.$region.hapcut2.threshold_30.output $chrom
/mendel-nas1/dhooper/SNPs/phased_vcfs/n50_extract.sh $tmp/$samid.$region.hapcut2.threshold_30.bed

##
echo "Preparing to concatenate phased variants and invariant site VCFs"
$bin_dir/bgzip $tmp/$samid.$region.hapcut2.threshold_30.output.phased.VCF
$bin_dir/tabix $tmp/$samid.$region.hapcut2.threshold_30.output.phased.VCF.gz
$bin_dir/bgzip $tmp/$sampleVCF
$bin_dir/tabix $tmp/$sampleVCF.gz

echo "Concatenating phased variants and invariant site VCFs"
$bin_dir/bcftools concat -a -d all $tmp/$samid.$region.hapcut2.threshold_30.output.phased.VCF.gz $tmp/$sampleVCF.gz -Oz -o $out_dir/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}
$bin_dir/tabix $out_dir/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}

echo "Cleaning up"
rm $tmp/$sampleHetVCF
rm $tmp/$sampleVCF.gz
rm $tmp/$samid.$region.unlinked.fragments
rm $tmp/$samid.$region.linked.fragments
rm $tmp/$samid.$region.hapcut2.threshold_30.output.phased.VCF.gz
rm $tmp/$samid*tbi

mv -v $tmp/$samid.$region.hapcut2.threshold_30.output $out_dir/
mv -v $tmp/$samid.$region.hapcut2.threshold_30.*bed $out_dir/