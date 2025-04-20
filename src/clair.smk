configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
sampleDir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing/2021.04.16'] # Check folder
refs = ['17762-33892_1_71_contigs.fa']
model = ['/home/belson/data/clair/ont/model']
rule all:
	input:
		expand('{results}/{sample}CorrectedPolca',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}.sorted.bam',results=results,sample=config["nanopore"]),
		#expand('{results}/{sample}CallVar',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}ClairPolca',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}_flye.clair.fasta',results=results,sample=config["nanopore"])
rule flye:
        input:
                expand('{sampleDir}/{sample}.fastq.gz',sampleDir=sampleDir,sample=config["nanopore"])
        output:
                directory('{results}/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids && samtools faidx {output}/assembly.fasta'
rule flyePolca:
	input:
		gen = rules.flye.output,
		r1 = expand('{sampleDir}/{r1}',sampleDir=sampleDir,r1=config["illumina"][0]),
                r2 = expand('{sampleDir}/{r2}',sampleDir=sampleDir,r2=config["illumina"][1])
	output:
		directory('{results}/{sample}CorrectedPolca')
	shell:
		"polca.sh -a {input.gen}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"
rule minimap:
	input:
		rules.flye.output
	output:
		temp('{results}/{sample}.sam')
	shell:
		'minimap2 -ax map-ont {input}/assembly.fasta {rules.flye.input} > {output}'

rule sortBam:
	input:
		rules.minimap.output
	output:
		temp('{results}/{sample}.sorted.bam')
	shell:
		'samtools sort -@ 8 -o {output} {input} && samtools index {output}'
rule clair:
	input:
		rules.sortBam.output
	output:
		'{results}/{sample}_flye.clair.fasta'
	shell:
		'./clair.sh {wildcards.sample} {wildcards.results}'
rule clairPolca:
	input:
		gen = rules.clair.output,
		r1 = expand('{sampleDir}/{r1}',sampleDir=sampleDir,r1=config["illumina"][0]),
                r2 = expand('{sampleDir}/{r2}',sampleDir=sampleDir,r2=config["illumina"][1])
	output:
		directory('{results}/{sample}ClairPolca')
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}_flye.clair.fasta* {output}"	
	
