configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
#samples = ['CHF10J']
root_dir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
res = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.03.24']
rule all:
	input:
		expand('{res}/{sample}_polish_flye_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_raconX1_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_raconX2_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_raconX3_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_raconX4_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_medaka_results',sample=config["nanopore"],res=res),
		expand('{res}/{sample}_polish_medaka2_results',sample=config["nanopore"],res=res),
		expand('{res}/polish_illumina_results',res=res)
		
rule denovo:
        input:
                expand('{root}/{sample}.fastq.gz',root=root_dir,sample=config["nanopore"])
        output:
                temp(directory('{res}/{sample}_flye'))
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polish_flye:
	input:
		gen = rules.denovo.output,
		r1 = expand('{root}/{read1}',root=root_dir,read1=config["illumina"][0]),
		r2 = expand('{root}/{read2}',root=root_dir,read2=config["illumina"][1])
	output:
		directory('{res}/{sample}_polish_flye_results')
	shell:
		"polca.sh -a {input.gen}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule raconX1:
	input:
		rules.denovo.output
	output:
		x1 = temp('{res}/{sample}_raconX1.fasta'),
		pf1 = temp('{res}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input}/assembly.fasta {rules.denovo.input} > {output.pf1} && racon -t 4 {rules.denovo.input} {output.pf1} {input}/assembly.fasta > {output.x1}'

rule polish_raconX1:
	input:
		rules.raconX1.output.x1
	output:
		directory('{res}/{sample}_polish_raconX1_results')
	shell:
		"polca.sh -a {input} -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv {wildcards.sample}_raconX1.fasta* {output}"

rule raconX2:
        input:
                rules.raconX1.output.x1
        output:
                x2 = temp('{res}/{sample}_raconX2.fasta'),
		pf2 = temp('{res}/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf2} && racon -t 4 {rules.denovo.input} {output.pf2} {input} > {output.x2}'

rule polish_raconX2:
	input:
		rules.raconX2.output.x2
	output:
		directory('{res}/{sample}_polish_raconX2_results')
	shell:
		"polca.sh -a {input} -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv {wildcards.sample}_raconX2.fasta* {output}"
rule raconX3:
        input:
                rules.raconX2.output.x2
        output:
                x3 = temp('{res}/{sample}_raconX3.fasta'),
		pf3 = temp('{res}/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf3} && racon -t 4 {rules.denovo.input} {output.pf3} {input} > {output.x3}'

rule polish_raconX3:
	input:
		rules.raconX3.output.x3
	output:
		directory('{res}/{sample}_polish_raconX3_results')
	shell:
		"polca.sh -a {input} -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv {wildcards.sample}_raconX3.fasta* {output}"
rule raconX4:
        input:
                rules.raconX3.output.x3
        output:
                x4 = ('{res}/{sample}_raconX4.fasta'),
		pf4 = temp('{res}/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.denovo.input} > {output.pf4} && racon -t 4 {rules.denovo.input} {output.pf4} {input} > {output.x4}'

rule polish_raconX4:
	input:
		rules.raconX4.output.x4
	output:
		directory('{res}/{sample}_polish_raconX4_results')
	shell:
		"polca.sh -a {input} -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv {wildcards.sample}_raconX4.fasta* {output}"
rule medaka:
	input:
		rules.raconX4.output.x4
	output:
		temp(directory('{res}/{sample}_medaka'))
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input} -t 8  -m r941_min_high_g303 -o {output}'

rule polish_medaka:
        input:
                rules.medaka.output
        output:
                directory('{res}/{sample}_polish_medaka_results')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
rule medaka2:
	input:
		rules.medaka.output
	output:
		temp(directory('{res}/{sample}_medaka2'))
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input}/consensus.fasta -t 8  -m r941_min_high_g303 -o {output}'

rule polish_medaka2:
        input:
                rules.medaka2.output
        output:
                directory('{res}/{sample}_polish_medaka2_results')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
rule spades:
	input:
		R1 = rules.polish_flye.input.r1,
		R2 = rules.polish_flye.input.r2
	output:
		temp(directory('{res}/illumina_results'))
	shell:
		'spades.py -1 {input.R1} -2 {input.R2} -o {output}'

rule polish_spades:
	input:
		rules.spades.output
	output:
		directory('{res}/polish_illumina_results')
	shell:
		"polca.sh -a {input}/contigs.fasta -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv contigs.fasta* {output}"
