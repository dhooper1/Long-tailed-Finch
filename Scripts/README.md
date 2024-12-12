# Scripts used for data analysis and their usage

## Linked-read mapping

Linked reads include BX tag information and require slight modification to default [BWA mem](https://github.com/lh3/bwa) mapping behavior. Also - see details [here](https://github.com/evolgenomics/haplotagging/tree/master) for recommended read mapping using barcode-first read mapper [EMA](https://github.com/arshajii/ema).

In the example below, the first 24 samples present in file (i.e., P04.sample.info) linking sample ID and well ID for a 96-well plate of samples are sequentially mapped to the [zebra finch reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003957565.4/).
```
# Map sequence data to zebra finch reference genome using BWA mem including BX tag information

reference="bTaeGut1.pri.cur.20210409.fasta"

for i in 01 02 03 04 05 06 07 08 09 `seq 10 24`; do
        name=$(grep C"$i" P04.sample.info | cut -f1);
        bwa mem -C -t 20 $reference \
        "$name"_C"$i"_P04_L003_NOH.R1.fastq.gz "$name"_C"$i"_P04_L003_NOH.R2.fastq.gz \
        -R "@RG\tID:"$name"\tSM:"$name"\tLB:"$name"\tPL:Illumina.NovaSeq4000.2x150" | samtools view -bh - > "$name"_C"$i"_P04_L003.bam;
        samtools sort -@ 20 -l 9 -T "$name".tmpsort -o "$name"_C"$i"_P04_L003.sorted.bam "$name"_C"$i"_P04_L003.bam;
        echo "Finished mapping and sorting sample: '$name'";
done
```

### Marking duplicates with BX-aware settings

```
# Mark duplicates in sorted BAM files with BX-aware options

for i in 01 02 03 04 05 06 07 08 09 `seq 10 24`; do
        name=$(grep C"$i" P04.sample.info | cut -f1);
        picard MarkDuplicates I="$name"_C"$i"_P04_L003.sorted.bam O="$name"_C"$i"_P04_L003.mkdup.bam M="$name"_C"$i".mkdup.metrics CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT;
        echo "Finished marking duplicates for sample: '$name'";
done
```

## Molecule summaries

Linked-read (LR) performance was evaluated by analyzing: number of unique (BX) beadTags, mean molecules per (BX) beadTag, mean read number per molecule, molecule N50 length, molecule mean length, and molecule maximum length. Note that each BAM file must be sorted by BX tag rather than by chromosome:position information.

Example approach for evaluating LR performance for sample AQ15:
```
# Sort a BAM file by BX tag information rather than by chrom:position
samtools sort -t BX AQ15_C02_P02.mkdup.bam -o AQ15_C02_P02.mkdup.BX-sorted.bam

# Write molecule output to a bed file
perl ./bed_write.pl AQ15_C02_P02.mkdup.BX-sorted.bam

# Filter results by mapping QUAL, beadTag specificity, and minimum number of reads per molecule
./filter_molecules.sh AQ15.P02.BX_sorted.linked_reads.full.bed

# Calculate molecule N50
./calculate_mol_N50.sh AQ15.P02.BX_sorted.filter.bed

```

## Variant calling

Initial variant calling performed using [bcftools](https://samtools.github.io/bcftools/bcftools.html) mpileup by chromosome using all project samples. Quality filtering of intial variant set performed using bcftools to generate a set of high-quality variants for imputation.

### Step #1: Use bcftools mpileup to create initial set of variants relative to the reference
In the example below, an intial set of SNP and INDEL variants are called relative to the zebra finch reference genome based on reads  from all project samples (i.e., *.mkdup.bam) mapped to a 5 Mb window on chromosome 2 (i.e., chr2).

```
reference="bTaeGut1.pri.cur.20210409.fasta"

bcftools mpileup -a INFO/AD -Ou -f $reference -r chr2:10000001-15000000 *.mkdup.bam | bcftools call -mv -f GQ,GP -Ob -o chr2:10000001-15000000.bcf
```

The variant calling process above was repeated in 5 Mb sliding-windows across the entire chromosome/genome.

### Step #2: Concatenate variant call 'chunks' into a single chromosome set
Variant calling was performed in 5 Mb sliding-windows across the entire chromosome, now we need to concatenate these into a single file for the entire chromosome.

```
# Create a list file containing the names of all window result files
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
bcftools filter -sLowQual -g3 -G10 -e'%QUAL<100 || (RPB<0.1 && %QUAL<50) || (AC<2 && %QUAL<50) || %MAX(AD)<=3 || %MAX(AD)/%MAX(DP)<=0.3' chr2.concatenated.bcf -Ob -o chr2.concatenated.filter.bcf

bcftools index chr2.concatenated.filter.bcf

bcftools view -m2 -M2 -v snps -i'FILTER="PASS"' chr2.concatenated.filter.bcf -Oz -o chr2.concatenated.biallelic.vcf.gz

bcftools index chr2.concatenated.biallelic.vcf.gz
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

```
# Create a list file of all STITCH output for a given chromosome
ls stitch.chr2.*vcf.gz > concat.chr2.list

# Concatenate results from a chromosome into a single file 
bcftools concat -f concat.chr2.list -a -D -Oz -o stitch.chr2.repeatmask.vcf.gz

bcftools index stitch.chr2.repeatmask.vcf.gz

# Filter variants with a STITCH INFO_SCORE <0.4
bcftools view -i'INFO_SCORE >= 0.4' stitch.chr2.repeatmask.vcf.gz -Oz -o stitch.chr2.repeatmask.filter.vcf.gz

bcftools index stitch.chr2.repeatmask.filter.vcf.gz

# Add the PG and PL fields to STITCH output
./add_PG_PL.sh stitch.chr2.repeatmask.filter.vcf.gz
```

## Imputation performance

Imputation performance was measured as the squared Pearson correlation between validation (i.e., high coverage) genotypes (GT) and imputed genotypes (GT) and dosages (DS) produced by STITCH. As imputation performance is strongly influenced by allele frequency, we measured it across eight representative frequency bins: [0-0.001], [0.001-0.002], [0.002-0.005], [0.005-0.01], [0.01-0.05], [0.05-0.1], [0.1-0.2], and [0.2-0.5]. Allele frequency was folded prior to evaluating imputation perforamnce so that no variant had an allele frequency that exceeded 0.5.

To evaluate imputation accuracy at the genotype level, genotype discordance was quantified between validation and imputed genotypes for: (1) all sites together, (2) homozygous major alleles, (3) heterozygous, and (4) homozygous minor alleles.

### Step #1: Prepare input data

In the example below, validation and imputed information is extracted for sample 'AG03' in order to evaluate imputation performance on chromosome 8.

```
# Extract validation genotype calls from desired sample of interest
bcftools query -s AG03 -f '%CHROM\t%POS\t[%GT]\n' chr8.validation.vcf.gz | gzip > AG03.chr8.validation.tsv.gz

# Extract STITCH INFO_SCORE, estimated allele frequency, genotype, and dosage for sample of interest
bcftools query -s AG03 -f '%CHROM\t%POS\t[%INFO_SCORE]\t[%EAF]\t[%GT]\t[%DS]\n' stitch.chr8.repeatmask.filter.vcf.gz | gzip > AG03.chr8.imputed.tsv.gz
```

### Step #2: Evaluate imputation performance

In the example below, imputation performance is evaluated as the Pearson correlation between validation and imputed genotypes (GT) and dosage (DS) at pre-defined minor allele frequency bins and as the genotype discordance at the genotype level. Results are written to screen and saved to desired output file.

```
# Imputation performance by genotype (GT) with a desired INFO_SCORE filter (i.e., 0.4)
python imputation_accuracy_GT_write2out.py AG03.chr8.validation.tsv.gz AG08.chr8.imputed.tsv.gz 0.4 AG03.chr8.GT.results.tsv

# Imputation performance by dosage (DS) with a desired INFO_SCORE filter (i.e., 0.4)
python imputation_accuracy_DS_write2out.py AG03.chr8.validation.tsv.gz AG08.chr8.imputed.tsv.gz 0.4 AG03.chr8.DS.results.tsv

# Imputation accuracy measured as genotype discordance with a desired INFO_SCORE and minor allele frequency filter
python genotype_discordance.py AG03.chr8.validation.tsv.gz AG08.chr8.imputed.tsv.gz 0.4 0.0
```

## Phasing using HapCUT2

Create a combined, phased VCF file for a single individual sample using: hapcutVcf.haplotagging.sh

In the example below, a text file (samples.list) providing sample IDs is used to phase all variants on chromosome 2 (i.e., chr2) given a multi-sample VCF file - generated above - containing each of the desired samples. Phasing carried out using [HapCUT2](https://github.com/vibansal/HapCUT2).

```
chrom="chr2"

cat samples.list | while read LINE; do
  ./hapcutVcf.haplotagging.sh ${LINE}_*.bam stitch.$chrom.repeatmask.filter.PL.vcf.gz $chrom;
done
```

### Calculate phase block N50 for a given chromosome

In the example below, phase block N50 is calculated on chr2 for all samples. 
```
chrom="chr2"
work_dir="/mendel-nas1/dhooper/SNPs/phased_vcfs/$chrom"
date="240329"

# Create a list of samples with phase block data
ls $work_dir/*n50.bed > $date.$chrom.n50.list

# Use calculate_phase_N50.sh script on each sample
cat $date.$chrom.n50.list | while read LINE; do
  ./calculate_phase_N50.sh ${LINE} >> $date.$chrom.n50.out;
done

rm $date.$chrom.n50.list
```

### Merge all phased single-sample VCF files into a multi-sample VCF

```
chrom="chr2"
date="231214"
work_dir="/mendel-nas1/dhooper/SNPs/phased_vcfs"

# Create a list of all VCF files you want to merge
ls *pMarkdup.$chrom.PL.AD.HAPCUT2.vcf.gz > $date.samples.$chrom.list

# Merge all single-sample VCF files in list file
bcftools merge -l $date.samples.$chrom.list -Oz -o hapcut2.$date.$chrom.merge.vcf.gz

tabix hapcut2.$date.$chrom.merge.vcf.gz
```

## Predict the landscape of recombination using ReLERNN

Ancestral recombination graph (ARG) analyses require a genomic map (i.e., the recombination landscape). In the example below, we demonstrate recombination landscape estimate using recurrent neural networks [ReLERNN](https://github.com/kr-colab/ReLERNN), a deep learning method for estimating a genome-wide recombination map using data from individually sequenced genomes.

### Step #1: Generate a VCF file for specific chromosome for a population of interest

In the example below we are focusing on chromosome 8 and only using sites with minor allele frequency >=2% and with a STITCH INFO_SCORE >=0.5. We are also only selecting a subset of male samples (N = 70) from subspecies *acuticauda* to match the sample set we would use to infer the recombination map of chromosome Z. 

```
out_dir="/mendel-nas1/dhooper/recombination/"
chr="chr8"

awk '{print $1}' LTFs.paa_allopatric_males.pop > tmp.list
bcftools view -S tmp.list stitch.$chr.repeatmask.filter.PL.vcf.gz | bcftools filter -i'MAF > 0.02 && INFO_SCORE>=0.5' -Oz -o $out_dir/paa_allo_males.$chr.repeatmask.filter.PL.MAF02.INFO_SCORE_0.5.vcf.gz
rm tmp.list

## ReLERNN requires input VCF be unzipped
cd $out_dir
gunzip paa_allo_males.$chr.repeatmask.filter.PL.MAF02.INFO_SCORE_0.5.vcf.gz
```

### Step #2: 

In the example below, we run through a full ReLERNN analysis to infer the recombination map on chromosome 8.

```
SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="5.85e-9"
URTR="15"
WORK_DIR="/mendel-nas1/dhooper/recombination"
DIR="$WORK_DIR/paa/chr8/"
VCF="$WORK_DIR/paa_allo_males.chr8.repeatmask.filter.PL.MAF02.INFO_SCORE_0.5.vcf"
GENOME="$WORK_DIR/genome.bed"
MASK="$WORK_DIR/chr8_accessibility_mask.bed"

##Unzip VCF file if currently zipped
gunzip $VCF.gz

module load Python/python-3.8.5-openmpi-cuda

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --mask ${MASK} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --assumedGenTime 1 \
    --upperRhoThetaRatio ${URTR} \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --nEpochs 500 \
    --nValSteps 500 \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --seed ${SEED}

conda activate base_genomics

bgzip $VCF

```

The file 'chr8.paa.map' was generated using results from the ReLERNN pipeline above.
