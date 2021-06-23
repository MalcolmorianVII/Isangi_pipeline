#! /usr/bin/bash

#NOTE:	$1 = {wildcards.sample}; $2 =  {wildcards.results}
set -e 
eval "$(conda shell.bash hook)"
conda activate "clair3"

model=/home/belson/Clair3/models/ont_guppy5
ref=$2/$1medaka/consensus.fasta
bam=$2/$1.sorted.bam

normVCF=$2/$1_flye_clair.norm.vcf.gz
samtools faidx $ref	#Somehow giving an error of not finding fai file for consensus.fasta
bash ~/Clair3/run_clair3.sh  --model_path=$model --ref_fn=$ref --bam_fn=$bam --platform="ont" --haploid_precise --output=$2 --threads=8 --include_all_ctgs
#Normalizing indels  & indexing the normalized vcf files
conda activate "bcf"
bcftools norm  -f $ref $2/merge_output.vcf.gz  -Oz -o $normVCF
tabix -p vcf $normVCF

#Altering the assembly
bcftools consensus -f $ref $normVCF > $2/$1_flye.clair.fasta
