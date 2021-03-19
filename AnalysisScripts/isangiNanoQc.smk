samples = ['FAO60983_pass_barcode02_b01ddb04']
refs = ['17762-33892_1_71']

rule all:
	input:
		expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_reads_stat.txt',sample=samples),
		expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_coverage.txt',sample=samples),
		expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_reads_mapped.txt',sample=samples),
		expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.11/{sample}_flye',sample=samples),
		expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.11/{sample}_assembly_stat.txt',sample=samples)
rule reads_stat:
	input:
		'/home/ubuntu/data/belson/test_data/2021.01.04/{sample}.fastq'
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_reads_stat.txt'
	shell:
		"assembly-stats -t {input} > {output}"
rule minimap:
	input:
		sam = rules.reads_stat.input,
		ref = expand('/home/ubuntu/data/belson/reference/2021.04.01/{ref}_contigs.fa',ref=refs)
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.sam')
	shell:
		'minimap2 -ax map-ont {input.ref} {input.sam} > {output}'

rule sam2bam:
	input:
		rules.minimap.output
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.bam')
	shell:
		'samtools view -b {input} -o {output}'

rule sort_bam:
	input:
		rules.sam2bam.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.sorted.bam'
	shell:
		'samtools sort {input} -o {output}'

rule sam_flagstat:
	input:
		rules.sort_bam.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_reads_mapped.txt'
	shell:
		'samtools flagstat {input} > {output}'

rule depth_calc:
	input:
		rules.sort_bam.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_coverage.txt'
	shell:
		"samtools depth -aa {input} | awk '{{sum+=$3}} END {{print \"Average = \",sum/NR}}' > {output}"

rule denovo:
	input:
		rules.reads_stat.input
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.11/{sample}_flye')
	shell:
		'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule assemb_stat:
	input:
		rules.denovo.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.11/{sample}_assembly_stat.txt'
	shell:
		"assembly-stats -t {input}/assembly.fasta > {output}"
