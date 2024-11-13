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

# Four rounds of genome polishing with Racon in a sequential order
rule raconX1:
	input:
		rules.denovo.output
	output:
		racon1 = '/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX1.fasta',
		paf1 = '/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon.paf'

	shell:
		'minimap2 -x map-ont {input}/assembly.fasta {rules.denovo.input} > {output.paf1} && racon -t 4 {rules.denovo.input} {output.paf1} {input}/assembly.fasta > {output.racon1}'
rule raconX2:
	# The input is the output of raconX1 stage
        input:
                rules.raconX1.output.racon1
        output:
                racon2 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX2.fasta'),
		paf2 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.paf2} && racon -t 4 {rules.denovo.input} {output.paf2} {input} > {output.racon2}'
rule raconX3:
	# The input is the output of raconX2 stage
        input:
                rules.raconX2.output.racon2
        output:
                racon3 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX3.fasta'),
		paf3 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.paf3} && racon -t 4 {rules.denovo.input} {output.paf3} {input} > {output.racon3}'
rule raconX4:
	# The input is the output of raconX3 stage
        input:
                rules.raconX3.output.racon3
        output:
                racon4 = ('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX4.fasta'),
		paf4 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.paf4} && racon -t 4 {rules.denovo.input} {output.paf4} {input} > {output.racon4}'
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
	
