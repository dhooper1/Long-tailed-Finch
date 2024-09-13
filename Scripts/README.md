### Scripts used for data analysis

## Variant calling

Initial variant calling performed using [bcftools](https://samtools.github.io/bcftools/bcftools.html) mpileup by chromosome using all project samples. Quality filtering of intial variant set performed using bcftools to generate a set of high-quality variants for imputation.

### Step #1: Use bcftools mpileup to create initial set of variants relative to the reference
In the example below, an intial set of SNP and INDEL variants are called relative to the [zebra finch reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003957565.4/) based on reads  from all project samples (i.e., *.mkdup.bam) mapped to a 5 Mb window on chromosome 2 (i.e., chr2).

```
reference="bTaeGut1.pri.cur.20210409.fasta"

bcftools mpileup -Ou -f $reference -r chr2:10000001-15000000 *.mkdup.bam | bcftools call -mv -Ob -o chr2:10000001-15000000.bcf
```

The variant calling process above was repeated in 5 Mb sliding-windows across the entire chromosome/genome.

### Step #2: Concatenate variant call 'chunks' into a single chromosome set
Variant calling was performed in 5 Mb sliding-windows across the entire chromosome, now we need to concatenate these into a single file for the entire chromosome.

```
ls chr2*.bcf > concat.chr2.list

# Index initial set of variants
cat concat.chr2.list | while read LINE; do
  bcftools index ${LINE};
done

# Concatenate all chunks into a single file
bcftools concat -f concat.chr2.list -a -D -Ob -o chr2.concatenated.bcf
bcftools index chr2.concatenated.bcf

```

### Step #3: Apply QUAL filters to initial variant call set and extract a set of biallelic SNP variants

```
bcftools filter -sLowQual -g3 -G10 -e'%QUAL<100 || (RPB<0.1 && %QUAL<50) || (AC<2 && %QUAL<50) || %MAX(AD)<=3 || %MAX(AD)/%MAX(DP)<=0.3' chr2.concatenated.bcf -Ob -o chr2.concatenated.filter

bcftools view -m2 -M2 -v snps -i'FILTER="PASS"' chr2.concatenated.filter.bcf -Oz -o chr2.concatenated.biallelic.vcf.gz
```

### Step #4: Remove SNP variants that overlap repeat masked regions
A BED file containing repeat regions called relative to the zebra finch reference genome is used to identify SNPs mappign to repetitive regions that we would like to remove.

```
repeats="bTaeGut1.pri.cur.20210409.fasta.repeatmask.bed"

vcftools --gzvcf chr2.concatenated.biallelic.vcf.gz --exclude-bed $repeats --recode --stdout | bgzip -c > chr2.concatenated.biallelic.repeatmask.vcf.gz

bcftools index chr2.concatenated.biallelic.repeatmask.vcf.gz
```

### Step #5: Extract SNP allelic details for STITCH

```
./make.position.file.STITCH.sh chr2
```

## Genotype imputation using STITCH

Imputation performed on a set of high-quality variants - generated above - using [STITCH](https://github.com/rwdavies/STITCH).

In the example below, characteristic for windows on chromosomes >25 Mb, imputation is performed on all desired variants in position file 'chr2.Q500.repeatmask.pos.txt' located within a 5 Mb window using reads from BAMs for all project samples listed in 'samples.bamlist.txt'.
```
R -e 'library("STITCH"); STITCH(tempdir = tempdir(), chr="chr2", bamlist="samples.bamlist.txt", posfile="chr2.Q500.repeatmask.pos.txt", outputdir="/mendel-nas1/dhooper/SNPs/STITCH/", method="pseudoHaploid", K=100, nGen=1000, S=1, readAware=TRUE, keepInterimFiles=FALSE, shuffle_bin_radius=100, iSizeUpperLimit=500000, keepSampleReadsInRAM=TRUE, niterations=40, switchModelIteration=25, regionStart=10000001, regionEnd=15000000, buffer=50000, expRate=1.0, outputSNPBlockSize=5000, use_bx_tag=FALSE, nCores=1)'
```

In the example, below, modifications have been made for 1 Mb windows on micro-chromosomes <10 Mb with much greater per-bp rates of recombination.
```
R -e 'library("STITCH"); STITCH(tempdir = tempdir(), chr="chr25", bamlist="samples.bamlist.txt", posfile="chr25.Q500.repeatmask.pos.txt", outputdir="/mendel-nas1/dhooper/SNPs/STITCH/", method="pseudoHaploid", K=100, nGen=1000, S=1, readAware=TRUE, keepInterimFiles=FALSE, shuffle_bin_radius=100, iSizeUpperLimit=500000, keepSampleReadsInRAM=TRUE, niterations=40, switchModelIteration=25, regionStart=1, regionEnd=1000000, buffer=50000, expRate=10.0, outputSNPBlockSize=5000, use_bx_tag=FALSE, nCores=1)'
```

The imputation approaches above were repeated in 5 Mb sliding-windows across the entire chromosome/genome before concatenating results for each chromosome.

## Phasing using HapCUT2

Create a combined, phased VCF file for a single individual sample using: hapcutVcf.haplotagging.sh

In the example below, a text file (samples.list) providing sample IID values is used to phase all variants on chromosome 2 (i.e., chr2) given a multi-sample VCF file containing each of the desired samples. Phasing carried out using [HapCUT2](https://github.com/vibansal/HapCUT2).

```
chrom="chr2"

cat samples.list | while read LINE; do
  ./hapcutVcf.haplo_230403.sh ${LINE}_*.bam stitch.$chrom.repeatmask.filter.PL.vcf.gz $chrom;
done
```
