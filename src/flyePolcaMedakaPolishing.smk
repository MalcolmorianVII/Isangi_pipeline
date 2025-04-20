configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
sampleDir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing/flye_polca_medaka']
rule all:
	input:
		expand('{results}/{sample}polishMedakaResults',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}polishFlyeResults',results=results,sample=config["nanopore"])

rule flye:
	input:
		expand('{sampleDir}/{sample}.fastq.gz',sampleDir=sampleDir,sample=config["nanopore"])
	output:
		directory('{results}/{sample}_flye')
	shell:
		'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polishFlye:
	input:
		assembly = rules.flye.output,
		r1 = expand('{sampleDir}/{r1}',sampleDir=sampleDir,r1=config["illumina"][0]),
		r2 = expand('{sampleDir}/{r2}',sampleDir=sampleDir,r2=config["illumina"][1])
	output:
		directory('{results}/{sample}polishFlyeResults')
	shell:
		"polca.sh -a {input.assembly}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule medaka:
        input:
                rules.polishFlye.output
        output:
               directory('{results}/{sample}Medaka')
        conda:
                '/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
        shell:
                'medaka_consensus -i {rules.flye.input} -d {input}/assembly.fasta.PolcaCorrected.fa -t 8  -m r941_min_high_g360 -o {output}'

rule polishMedaka:
        input:
                rules.medaka.output
        output:
                directory('{results}/{sample}polishMedakaResults')
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polishFlye.input.r1} {rules.polishFlye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
