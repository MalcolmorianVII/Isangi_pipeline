configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
sampleDir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.03.24']
rule all:
	input:
		expand('{results}/{sample}polishFlye',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishRaconX1',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishRaconX2',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishRaconX3',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishRaconX4',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishMedaka',sample = config["nanopore"],results = results),
		expand('{results}/{sample}polishMedaka2',sample = config["nanopore"],results = results),
		expand('{results}/polishIllumina',results = results)
		
rule flye:
        input:
                expand('{sampleDir}/{sample}.fastq.gz',sampleDir = sampleDir,sample = config["nanopore"])
        output:
                temp(directory('{results}/{sample}Flye'))
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polishFlye:
	input:
		gen = rules.flye.output,
		r1 = expand('{sampleDir}/{read1}',sampleDir = sampleDir,read1=config["illumina"][0]),
		r2 = expand('{sampleDir}/{read2}',sampleDir = sampleDir,read2=config["illumina"][1])
	output:
		directory('{results}/{sample}polishFlye')
	shell:
		"polca.sh -a {input.gen}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule raconX1:
	input:
		rules.flye.output
	output:
		x1 = '{results}/{sample}RaconX1.fasta',
		pf1 = temp('{results}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input}/assembly.fasta {rules.flye.input} > {output.pf1} && racon -t 4 {rules.flye.input} {output.pf1} {input}/assembly.fasta > {output.x1}'

rule polish_raconX1:
	input:
		rules.raconX1.output.x1
	output:
		directory('{results}/{sample}polishRaconX1')
	shell:
		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX1.fasta* {output}"

rule raconX2:
        input:
                rules.raconX1.output.x1
        output:
                x2 = temp('{results}/{sample}RaconX2.fasta'),
		pf2 = temp('{results}/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.flye.input} > {output.pf2} && racon -t 4 {rules.flye.input} {output.pf2} {input} > {output.x2}'

rule polish_raconX2:
	input:
		rules.raconX2.output.x2
	output:
		directory('{results}/{sample}polishRaconX2')
	shell:
		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX2.fasta* {output}"
rule raconX3:
        input:
                rules.raconX2.output.x2
        output:
                x3 = temp('{results}/{sample}RaconX3.fasta'),
		pf3 = temp('{results}/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.flye.input} > {output.pf3} && racon -t 4 {rules.flye.input} {output.pf3} {input} > {output.x3}'

rule polish_raconX3:
	input:
		rules.raconX3.output.x3
	output:
		directory('{results}/{sample}polishRaconX3')
	shell:
		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX3.fasta* {output}"
rule raconX4:
        input:
                rules.raconX3.output.x3
        output:
                x4 = ('{results}/{sample}RaconX4.fasta'),
		pf4 = temp('{results}/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} {rules.flye.input} > {output.pf4} && racon -t 4 {rules.flye.input} {output.pf4} {input} > {output.x4}'

rule polish_raconX4:
	input:
		rules.raconX4.output.x4
	output:
		directory('{results}/{sample}polishRaconX4')
	shell:
		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX4.fasta* {output}"
rule medaka:
	input:
		rules.raconX4.output.x4
	output:
		temp(directory('{results}/{sample}Medaka'))
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.flye.input} -d {input} -t 8  -m r941_min_high_g303 -o {output}'

rule polish_medaka:
        input:
                rules.medaka.output
        output:
                directory('{results}/{sample}polishMedaka')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
rule medaka2:
	input:
		rules.medaka.output
	output:
		temp(directory('{results}/{sample}Medaka2'))
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {rules.flye.input} -d {input}/consensus.fasta -t 8  -m r941_min_high_g303 -o {output}'

rule polish_medaka2:
        input:
                rules.medaka2.output
        output:
                directory('{results}/{sample}polishMedaka2')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
rule spades:
	input:
		R1 = rules.polishFlye.input.r1,
		R2 = rules.polishFlye.input.r2
	output:
		temp(directory('{results}/illuminaResults'))
	shell:
		'spades.py -1 {input.R1} -2 {input.R2} -o {output}'

rule polish_spades:
	input:
		rules.spades.output
	output:
		directory('{results}/polishIllumina')
	shell:
		"polca.sh -a {input}/contigs.fasta -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv contigs.fasta* {output}"
