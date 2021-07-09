configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
samples = ['1_Acinetobacter_baumannii_J9','2_Citrobacter_koseri_MINF_9D','3_Enterobacter_kobei_MSB1_1B','4_Haemophilus_unknown_M1C132_1','5_Klebsiella_oxytoca_MSB1_2C','7_Klebsiella_variicola_INF345','8_Serratia_marcescens_17-147-1671']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing/2021.07.08/Other_species/Guppy3']
nano = '/home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads'
#def lambda wildcards : config[wildcards.sample]["R1"](wildcards):
#	return expand('{read1}',read1=config[wildcards.sample]["R1"])
#def  lambda wildcards : config[wildcards.sample]["R2"](wildcards):
#	return expand('{read2}',read2=config[wildcards.sample]["R2"])
rule all:
	input:
		expand('{results}/{sample}/{sample}Flye',sample = samples,results = results),
		expand('{results}/{sample}/{sample}FlyePolished',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishRaconX1',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishRaconX2',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishRaconX3',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishRaconX4',sample = samples,results = results),
		expand('{results}/{sample}/{sample}medaka',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishMedaka',sample = samples,results = results),
		#expand('{results}/{sample}/{sample}polishMedaka2',sample = samples,results = results),
		expand('{results}/{sample}/{sample}polishIllumina',results = results,sample=samples)
		
rule flye:
	input:
		'/home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{sample}.fastq.gz'	# {{sample}} 2 indicate it is a wildcard from output
	output:
                directory('{results}/{sample}/{sample}Flye')
	shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polishFlye:
	input:
		gen = rules.flye.output,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
		r2 = lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}FlyePolished')
	shell:
		"polca.sh -a {input.gen}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule raconX1:
	input:
		rules.flye.output
	output:
		x1 = temp('{results}/{sample}/{sample}RaconX1.fasta'),
		pf1 = temp('{results}/{sample}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input}/assembly.fasta /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz > {output.pf1} && racon -t 4 /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz {output.pf1} {input}/assembly.fasta > {output.x1}'

rule polish_raconX1:
	input:
		gen = rules.raconX1.output.x1,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
		r2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}polishRaconX1')
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX1.fasta* {output}"

rule raconX2:
        input:
                rules.raconX1.output.x1
        output:
                x2 = temp('{results}/{sample}/{sample}RaconX2.fasta'),
		pf2 = temp('{results}/{sample}/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input} /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz > {output.pf2} && racon -t 4 /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz {output.pf2} {input} > {output.x2}'

rule polish_raconX2:
	input:
		gen = rules.raconX2.output.x2,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
                r2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}polishRaconX2')
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX2.fasta* {output}"
rule raconX3:
        input:
                rules.raconX2.output.x2
        output:
                x3 = temp('{results}/{sample}/{sample}RaconX3.fasta'),
		pf3 = temp('{results}/{sample}/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input} /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz > {output.pf3} && racon -t 4 /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz {output.pf3} {input} > {output.x3}'

rule polish_raconX3:
	input:
		gen = rules.raconX3.output.x3,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
                r2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}polishRaconX3')
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX3.fasta* {output}"
rule raconX4:
        input:
                rules.raconX3.output.x3
        output:
                x4 = ('{results}/{sample}/{sample}RaconX4.fasta'),
		pf4 = temp('{results}/{sample}/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input} /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz > {output.pf4} && racon -t 4 /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz {output.pf4} {input} > {output.x4}'

rule polish_raconX4:
	input:
		gen = rules.raconX4.output.x4,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
                r2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}polishRaconX4')
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}RaconX4.fasta* {output}"
rule medaka:
	input:
		rules.raconX4.output.x4
	output:
		directory('{results}/{sample}/{sample}medaka')
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i /home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads/barcode0{wildcards.sample}.fastq.gz -d {input} -t 8  -m r941_min_high_g360 -o {output}'

rule polish_medaka:
        input:
                gen = rules.medaka.output,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
                r2 =  lambda wildcards : config[wildcards.sample]["R2"]
        output:
                directory('{results}/{sample}/{sample}polishMedaka')
        shell:
                "polca.sh -a {input.gen}/consensus.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv consensus.fasta* {output}"

rule spades:
	input:
		R1 = lambda wildcards : config[wildcards.sample]["R1"],
		R2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}illuminaResults')
	shell:
		'spades.py -1 {input.R1} -2 {input.R2} --phred-offset 33 -o {output}'

rule polish_spades:
	input:
		gen = rules.spades.output,
		r1 = lambda wildcards : config[wildcards.sample]["R1"],
                r2 =  lambda wildcards : config[wildcards.sample]["R2"]
	output:
		directory('{results}/{sample}/{sample}polishIllumina')
	shell:
		"polca.sh -a {input.gen}/contigs.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv contigs.fasta* {output}"
