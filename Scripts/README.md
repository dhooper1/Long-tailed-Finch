### Scripts used for data analysis

## Phasing using HapCUT2

Create a combined, phased VCF file for a single individual sample using: hapcutVcf.haplotagging.sh

In the example below, a text file (samples.list) providing sample IID values is used to phase all variants on chromosome 2 (i.e., chr2) given a multi-sample VCF file containing each of the desired samples.

```
chrom="chr2"
cat samples.list | while read LINE; do
  ./hapcutVcf.haplo_230403.sh ${LINE}_*.bam stitch.$chrom.repeatmask.filter.PL.vcf.gz $chrom;
done
```
