configfile:"config.yaml"
#samples = ['CHF10J']
root_dir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
rule all:
        input:
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06/{sample}_circularised',sample=config["nanopore"]),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06/{sample}_amrfinder',sample=config["nanopore"]),
                expand('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06/{sample}_prokka',sample=config["nanopore"])

rule denovo:
        input:
                expand('{root}/{sample}.fastq.gz',root=root_dir,sample=config["nanopore"])
        output:
                directory('{root}/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polish_flye:
	input:
		rules.denovo.output
	output:
		directory('{root}/polish_flye_results')
	shell:
		"polca.sh -a {input}/assembly.fasta -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"

rule raconX1:
	input:
		rules.denovo.output
	output:
		racon1 = temp('{root}/{sample}_raconX1.fasta'),
		paf1 = temp('{root}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input}/assembly.fasta {rules.denovo.input} > {output.pf1} && racon -t 4 {rules.denovo.input} {output.pf1} {input}/assembly.fasta > {output.x1}'

rule polish_raconX1:
	input:
		rules.raconX1.output.racon1
	output:
		directory('{root}/polish_raconX1_results')
	shell:
		"polca.sh -a {input} -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"

rule raconX2:
        input:
                rules.raconX1.output.racon1
        output:
                racon2 = temp('{root}/{sample}_raconX2.fasta'),
		paf2 = temp('{root}/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf2} && racon -t 4 {rules.denovo.input} {output.pf2} {input} > {output.x2}'

rule polish_raconX2:
	input:
		rules.raconX2.output.racon2
	output:
		directory('{root}/polish_raconX2_results')
	shell:
		"polca.sh -a {input} -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"
rule raconX3:
        input:
                rules.raconX2.output.racon2
        output:
                racon3 = temp('{root}/{sample}_raconX3.fasta'),
		paf3 = temp('{root}/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf3} && racon -t 4 {rules.denovo.input} {output.pf3} {input} > {output.x3}'

rule polish_raconX3:
	input:
		rules.raconX3.output.racon3
	output:
		directory('{root}/polish_raconX3_results')
	shell:
		"polca.sh -a {input} -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"
rule raconX4:
        input:
                rules.raconX3.output.racon3
        output:
                racon4 = ('{root}/{sample}_raconX4.fasta'),
		paf4 = temp('{root}/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf4} && racon -t 4 {rules.denovo.input} {output.pf4} {input} > {output.x4}'

rule polish_raconX4:
	input:
		rules.raconX4.output.racon4
	output:
		directory('{root}/polish_raconX4_results')
	shell:
		"polca.sh -a {input} -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"
rule medaka:
	input:
		rules.raconX4.output.racon4
	output:
		directory('{root}/{sample}_medaka')
	conda:
		'envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input} -t 8  -m r941_min_fast_g303 -o {output}'

rule polish_medaka:
	input:
		rules.medaka.output
	output:
		directory('{root}/polish_medaka_results')
	shell:
		"polca.sh -a {input}/consensus.fasta -r'{wildcards.root}/17762-33892_1_71_bbduk_1.fastq.gz {wildcards.root}/17762-33892_1_71_bbduk_2.fastq.gz'"
rule circlator:
	input:
		rules.medaka.output
	output:
		directory('{root}/{sample}_circularised')
	shell:
		'circlator all --merge_min_id 85 --merge_breaklen 1000 {input}/*.fasta {rules.denovo.input} {output}'


rule amrfinder:
	input:
		rules.circlator.output

	output:
		'{root}/{sample}_amrfinder'
	shell:
		'amrfinder --plus -n {input}/06.fixstart.fasta -O Salmonella > {output}'

rule prokka:
	input:
		rules.circlator.output
	output:
		directory('{root}/{sample}_prokka')
	conda:
		'envs/prokka.yml'
	shell:
		'prokka {input}/06.fixstart.fasta --outdir {output}'
	
