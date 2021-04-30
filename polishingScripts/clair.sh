#!/usr/bin/bash

#NOTE:	$1 = {wildcards.sample}; $2 =  {wildcards.results}
set -e 
eval "$(conda shell.bash hook)"
conda activate "clair-env"

model=/home/belson/data/clair/ont/model
ref=$2/$1_flye/assembly.fasta
bam=$2/$1.sorted.bam
callVar=$2/$1CallVar
flyeClair=$2/$1_flye_clair.vcf
flyeClairgz=$flyeClair.gz
normVCF=$2/$1_flye_clair.norm.vcf.gz
clair.py callVarBamParallel --chkpnt_fn $model --ref_fn $ref --bam_fn $bam --threshold 0.2 --sampleName $1 --haploid_precision --includingAllContigs --output_prefix $callVar > $callVar.command.sh

cat $callVar.command.sh | parallel -j2

# prep vcf files for bcftools
vcfcat $2/*.vcf > $flyeClair
bgzip $flyeClair
tabix -p vcf $flyeClairgz

#Normalizing indels  & indexing the normalized vcf files
bcftools norm  -f $ref $flyeClairgz -m +any -Oz -o $normVCF
tabix -p vcf $normVCF

#Altering the assembly
bcftools consensus -f $ref $normVCF -H R > $2/$1_flye.clair.fasta
