samples = ['FAO60983_pass_barcode02_b01ddb04']
root_dir = ['/home/ubuntu/data/belson/test_data/2021.01.04']
rule all:
        input:
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_circularised',sample=samples),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder',sample=samples),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka',sample=samples)


rule denovo:
        input:
                expand('{root}/{sample}.fastq',root=root_dir,sample=samples)
        output:
                directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule raconX1:
	input:
		rules.denovo.output
	output:
		x1 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX1.fasta'),
		pf1 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input} {rules.denovo.input} > {ouput.pf1} && racon -t 4 {rules.denovo.input} {output.pf1} {input}/assembly.fasta > {output}'
rule raconX2:
        input:
                rules.raconX1.output
        output:
                x2 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX2.fasta'),
		pf2 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf2} && racon -t 4 {rules.denovo.input} {output.pf2} {input} > {output}'
rule raconX3:
        input:
                rules.raconX2.output
        output:
                x3 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX3.fasta'),
		pf3 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf3} && racon -t 4 {rules.denovo.input} {output.pf3} {input} > {output}'
rule raconX4:
        input:
                rules.raconX3.output
        output:
                x4 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_raconX4.fasta'),
		pf4 = temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf4} && racon -t 4 {rules.denovo.input} {output.pf4} {input} > {output}'
rule medaka:
	input:
		rules.raconX4.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_medaka')
	conda:
		'envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input} -t 8  -m r941_min_fast_g303 -o {output}'

rule circlator:
	input:
		rules.medaka.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_circularised')
	shell:
		'circlator all --merge_min_id 85 --merge_breaklen 1000 {input}/*.fasta {rules.denovo.input} {output}'


rule amrfinder:
	input:
		rules.circlator.output

	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_amrfinder'
	shell:
		'amrfinder --plus -n {input}/06.fixstart.fasta -O Salmonella > {output}'

rule prokka:
	input:
		rules.circlator.output
	output:
		directory('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.18/{sample}_prokka')
	conda:
		'envs/prokka.yml'
	shell:
		'prokka {input}/06.fixstart.fasta --outdir {output}'
	
