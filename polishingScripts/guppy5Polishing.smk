configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
samples = ['1_Acinetobacter_baumannii_J9','2_Citrobacter_koseri_MINF_9D','3_Enterobacter_kobei_MSB1_1B','4_Haemophilus_unknown_M1C132_1','5_Klebsiella_oxytoca_MSB1_2C','7_Klebsiella_variicola_INF345','8_Serratia_marcescens_17-147-1671']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing/2021.07.14/Other_species/Guppy5']
#nano = '/home/ubuntu/data/belson/monash_data/guppy_v5.0.7_reads'
root_dir = config['torun']
model = config['model']
def reads(wildcards):
	return expand('{root_dir}/barcode0{{sample}}.fastq.gz',root_dir=root_dir,sample=samples)
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
		nano=reads
	output:
		directory('{results}/{sample}/{sample}Flye')
	shell:
		'flye --nano-raw {input.nano} -g 5m -o {output} -t 8 --plasmids'

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
		gen = rules.flye.output,
		nano = reads
	output:
		x1 = temp('{results}/{sample}/{sample}RaconX1.fasta'),
		pf1 = temp('{results}/{sample}/{sample}.racon.paf')
	shell:
		'minimap2 -x map-ont {input.gen}/assembly.fasta {input.nano} > {output.pf1} && racon -t 4 {input.nano} {output.pf1} {input.gen}/assembly.fasta > {output.x1}'

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
                gen = rules.raconX1.output.x1,
		nano = reads
        output:
                x2 = temp('{results}/{sample}/{sample}RaconX2.fasta'),
		pf2 = temp('{results}/{sample}/{sample}.racon2.paf')
        shell:
                'minimap2 -x map-ont {input.gen} {input.nano} > {output.pf2} && racon -t 4 {input.nano} {output.pf2} {input.gen} > {output.x2}'

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
                gen = rules.raconX2.output.x2,
		nano = reads
        output:
                x3 = temp('{results}/{sample}/{sample}RaconX3.fasta'),
		pf3 = temp('{results}/{sample}/{sample}.racon3.paf')
        shell:
                'minimap2 -x map-ont {input.gen} {input.nano} > {output.pf3} && racon -t 4 {input.nano} {output.pf3} {input.gen} > {output.x3}'

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
                gen = rules.raconX3.output.x3,
		nano = reads
        output:
                x4 = ('{results}/{sample}/{sample}RaconX4.fasta'),
		pf4 = temp('{results}/{sample}/{sample}.racon4.paf')
        shell:
                'minimap2 -x map-ont {input.gen} {input.nano} > {output.pf4} && racon -t 4 {input.nano} {output.pf4} {input.gen} > {output.x4}'

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
		gen = rules.raconX4.output.x4,
		nano = reads
	output:
		directory('{results}/{sample}/{sample}medaka')
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
	shell:
		'medaka_consensus -i {input.nano} -d {input.gen} -t 8  -m {model} -o {output}'

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
