samples = ['FAO60983_pass_barcode02_b01ddb04']
rule all:
        input:
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_circularised',sample=samples),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder',sample=samples),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka',sample=samples)


rule denovo:
        input:
                '/home/ubuntu/data/belson/test_data/2021.01.04/{sample}.fastq'
        output:
                directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule raconX1:
	input:
		rules.denovo.output
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX1.fasta')
	shell:
		'racon {input}/assembly_fasta > {output}'
rule raconX2:
        input:
                rules.raconX1.output
        output:
                temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX2.fasta')
        shell:
                'racon {input} > {output}'
rule raconX3:
        input:
                rules.raconX2.output
        output:
                temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX3.fasta')
        shell:
                'racon {input} > {output}'
rule raconX4:
        input:
                rules.raconX3.output
        output:
                temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX4.fasta')
        shell:
                'racon {input} > {output}'
rule medaka:
	input:
		rules.raconX4.output
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_medaka')
	shell:
		'medaka consensus -d {input} -o {output}'

rule circlator:
	input:
		rules.medaka.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_circularised')
	shell:
		'circlator all --merge_min_id 85 --merge_breaklen 1000 {input}/*.fasta {output}'


rule amrfinder:
	input:
		rules.circlator.output

	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder')
	shell:
		'amrfinder --plus -n {input}/06.fixstart.fasta -O Salmonella > {output}'

rule prokka:
	input:
		rules.circlator.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka')
	shell:
		'prokka {input}/06.fixstart.fasta --outdir {output}'
	
