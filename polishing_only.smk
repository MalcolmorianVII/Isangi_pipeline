samples = ['CHF10J']
root_dir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
res = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.16']
rule all:
	input:
		#expand('{res}/{sample}_medaka',sample=samples,res=res),
		expand('{res}/{sample}_polish_flye_results',sample=samples,res=res),
		expand('{res}/{sample}_polish_raconX1_results',sample=samples,res=res),
		expand('{res}/{sample}_polish_raconX2_results',sample=samples,res=res),
		expand('{res}/{sample}_polish_raconX3_results',sample=samples,res=res),
		expand('{res}/{sample}_polish_raconX4_results',sample=samples,res=res),
		expand('{res}/{sample}_polish_medaka_results',sample=samples,res=res)

rule denovo:
        input:
                expand('{root}/{sample}.fastq.gz',root=root_dir,sample=samples)
        output:
                directory('{res}/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polish_flye:
	input:
		gen = rules.denovo.output,
		r1 = expand('{root}/17762-33892_1_71_bbduk_1.fastq.gz',root=root_dir),
		r2 = expand('{root}/17762-33892_1_71_bbduk_2.fastq.gz',root=root_dir)
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
		directory('{res}/{sample}_medaka')
	conda:
		'envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.denovo.input} -d {input} -t 8  -m r941_min_fast_g303 -o {output}'

rule polish_medaka:
        input:
                rules.medaka.output
        output:
                directory('{res}/{sample}_polish_medaka_results')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
