# Genomic analysis workflow i.e denovo assembly + genome polishing + annotation
samples = ['FAO60983_pass_barcode02_b01ddb04']
root_dir = ['/home/ubuntu/data/belson/test_data/2021.01.04']
rule all:
        input:
		# 2 results of interest i.e The AMRfinder results file & prokka annnotations
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder',sample=samples),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka',sample=samples)


rule denovo:
	# First assembling the genome i.e using flye
        input:
                expand('{root}/{sample}.fastq',root=root_dir,sample=samples)
        output:
                directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'


# Map the accurate short reads to the assembly for polishing
rule  minimap:
	input:
		rules.denovo.output,
		r1 = '/home/ubuntu/samples/{sample}_1.fastq',
		r2 = '/home/ubuntu/samples/{sample}_2.fastq'
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon.paf')

	shell:
		'minimap2 -x sr {input}/assembly.fasta {input.r1} {input.r2} > {output}'

# Genome polishing with Racon
rule racon:
	input:
		rules.denovo.output
	output:
		racon1 = '/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX1.fasta',
		paf1 = '/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon.paf'

	shell:
		'minimap2 -x map-ont {input}/assembly.fasta {rules.denovo.input} > {output.paf1} && racon -t 4 {rules.denovo.input} {output.paf1} {input}/assembly.fasta > {output.racon1}'

rule medaka:
	# The input is the output of raconX4 stage
	input:
		rules.raconX4.output.racon4
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_medaka')
	conda:
		'envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input} -t 8  -m r941_min_fast_g303 -o {output}'

rule circlator:
	# Circularizing the draft genome after polishing i.e coorecting breaks introduced during polishing & assmebly????
	input:
		rules.medaka.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_circularised')
	shell:
		'circlator all --merge_min_id 85 --merge_breaklen 1000 {input}/*.fasta {rules.denovo.input} {output}'


rule amrfinder:
	# Probing the AMRgenes in the draft genome
	input:
		rules.circlator.output

	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder'
	shell:
		'amrfinder --plus -n {input}/06.fixstart.fasta -O Salmonella > {output}'

rule prokka:
	# Annotating the genome
	input:
		rules.circlator.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka')
	conda:
		'envs/prokka.yml'
	shell:
		'prokka {input}/06.fixstart.fasta --outdir {output}'
	
