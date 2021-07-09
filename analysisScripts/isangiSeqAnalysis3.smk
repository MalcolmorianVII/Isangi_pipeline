configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/try.yaml"
#sampleDir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.05.15']
rule all:
        input:
                expand('{results}/{sample}Amrfinder',results = results,sample=config["samples"]),
                expand('{results}/{sample}Bakta',results = results,sample=config["samples"]),
#		expand('{results}/{sample}polishFlyeResults',results = results,sample=config["samples"])

rule flye:
        input:
                lambda wildcards: expand('{sample}',sample =config["samples"][wildcards.sample]["nanopore"])
        output:
                directory('{results}/{sample}Flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polishFlye:
	input:
                r1 = lambda wildcards: config["samples"][wildcards.sample]["illumina"][0],
                r2 = lambda wildcards: config["samples"][wildcards.sample]["illumina"][1],
		gen = rules.flye.output
	output:
		directory('{results}/{sample}polishFlyeResults')
	shell:
		"polca.sh -a {input.gen}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule raconX1:
	input:
		gen = rules.flye.output,
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
	output:
		racon1 = temp('{results}/{sample}RaconX1.fasta'),
		paf1 = temp('{results}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input.gen}/assembly.fasta {input.nano} > {output.paf1} && racon -t 4 {input.nano} {output.paf1} {input.gen}/assembly.fasta > {output.racon1}'

#rule polishRaconX1:
#	input:
#		rules.raconX1.output.racon1
#	output:
#		directory('{results}/polishRaconX1Results')
#	shell:
#		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX1.fasta* {output}"

rule raconX2:
	input:
		geno = rules.raconX1.output.racon1,
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
	output:
                racon2 = str(temp('{results}/{sample}RaconX2.fasta')),
		paf2 = str(temp('{results}/{sample}.racon2.paf'))
	shell:
                'minimap2 -x map-ont {input.geno} {input.nano} > {output.paf2} && racon -t 4 {input.nano} {output.paf2} {input.geno} > {output.racon2}'

#rule polishRaconX2:
#	input:
#		rules.raconX2.output.racon2
#	output:
#		directory('{results}/polishRaconX2Results')
#	shell:
#		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX2.fasta* {output}"
rule raconX3:
	input:
		geno = rules.raconX2.output.racon2,
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
	output:
                racon3 = str(temp('{results}/{sample}RaconX3.fasta')),
		paf3 = str(temp('{results}/{sample}.racon3.paf'))
	shell:
                'minimap2 -x map-ont {input.geno} {input.nano} > {output.paf3} && racon -t 4 {input.nano} {output.paf3} {input.geno} > {output.racon3}'

#rule polishRaconX3:
#	input:
#		rules.raconX3.output.racon3
#	output:
#		directory('{results}/polishRaconX3Results')
#	shell:
#		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX3.fasta* {output}"
rule raconX4:
	input:
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"],
		geno = rules.raconX3.output.racon3
	output:
                racon4 = str('{results}/{sample}RaconX4.fasta'),
		paf4 = str(temp('{results}/{sample}.racon4.paf'))
	shell:
                'minimap2 -x map-ont {input.geno} {input.nano} > {output.paf4} && racon -t 4 {input.nano} {output.paf4} {input.geno} > {output.racon4}'

#rule polishRaconX4:
#	input:
#		rules.raconX4.output.racon4
#	output:
#		directory('{results}/polishRaconX4Results')
#	shell:
#		"polca.sh -a {input} -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX4.fasta* {output}"
rule medaka:
	input:
		geno = rules.raconX4.output.racon4,
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
	output:
		directory('{results}/{sample}Medaka')
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {input.nano} -d {input.geno} -t 8  -m r941_min_high_g360 -o {output}'

#rule polishMedaka:
#	input:
#		rules.medaka.output
#	output:
#		directory('{results}/polishMedakaResults')
#	shell:
#		"polca.sh -a {input}/consensus.fasta -r'{rules.polishFlye.r1} {rules.polishFlye.r2}' && mkdir {output} && mv consensus.fasta* {output}"
rule circlator:
	input:
		geno = rules.medaka.output,
		nano = lambda wildcards: config["samples"][wildcards.sample]["nanopore"]
	output:
		directory('{results}/{sample}Circularised')
	shell:
		'circlator all --merge_min_id 85 --merge_breaklen 1000 {input.geno}/*.fasta {input.nano} {output}'


rule amrfinder:
	input:
		rules.circlator.output

	output:
		'{results}/{sample}Amrfinder'
	shell:
		'amrfinder --plus -n {input}/06.fixstart.fasta -O Salmonella > {output}'
rule bakta:
        input:
                rules.circlator.output
        output:
                directory('{results}/{sample}Bakta')
        shell:
                './bakta.sh {input}/06.fixstart.fasta {output}'
#rule prokka:
#	input:
#		rules.circlator.output
#	output:
#		directory('{results}/{sample}Prokka')
#	conda:
#		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/prokka.yml'
#	shell:
#		'prokka {input}/06.fixstart.fasta --prefix {wildcards.sample} --genus Salmonella --species enterica --outdir {output}'
	
