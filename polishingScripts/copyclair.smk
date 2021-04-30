configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
sampleDir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.04.16'] # Check folder
refs = ['17762-33892_1_71_contigs.fa']
model = ['/home/belson/data/clair/ont/model']
rule all:
	input:
		expand('{results}/{sample}.sorted.bam',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}CallVar',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}ClairPolca',results=results,sample=config["nanopore"]),
		expand('{results}/{sample}_flye.clair.fasta',results=results,sample=config["nanopore"])
rule flye:
        input:
                expand('{sampleDir}/{sample}.fastq.gz',sampleDir=sampleDir,sample=config["nanopore"])
        output:
                directory('{results}/{sample}_flye')
        shell:
                'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'
		
rule minimap:
	input:
		rules.flye.output
	output:
		'{results}/{sample}.sam'
	shell:
		'minimap2 -ax map-ont {input}/assembly.fasta {rules.flye.input} > {output}'

rule sortBam:
	input:
		rules.minimap.output
	output:
		'{results}/{sample}.sorted.bam'
	shell:
		'samtools sort -@ 8 -o {output} {input} && samtools index {output}'

rule callVarBamParallel:
	input:
		bam = rules.sortBam.output,
		ref = rules.flye.output
		#callVar = rules.CVMprep.output,
	output:
		callVar = directory('{results}/{sample}CallVar'),
		command = '{results}/{sample}CallVar/command.sh'
	conda:
		#directory('/home/belson/anaconda3/envs/clair-env/bin/clair')
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/clair.yml'
	shell:
		'clair.py callVarBamParallel --chkpnt_fn {model} --ref_fn {input.ref}/assembly.fasta --bam_fn {input.bam} --threshold 0.2 --sampleName {wildcards.sample} --haploid_precision --includingAllContigs --output_prefix {output.callVar} > {output.command}'
rule parallelize:
	input:
		rules.callVarBamParallel.output.command
	output:
		ct1 = '{results}/{sample}/CallVar.contig_1_0_4729009.vcf',
		ct2 = '{results}/{sample}/CallVar.contig_2_0_203003.vcf'
	conda:
		#directory('/home/belson/anaconda3/envs/clair-env/bin/clair')
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/clair.yml'
	shell:
		'cat {input} | parallel -j2'

rule vcfcat:
	input:
		ct1 = rules.parallelize.output.ct1,
		ct2 = rules.parallelize.output.ct2,
		#rules.callVarBamParallel.output.callVar
	output:
		'{results}/{sample}_flye_clair.vcf'
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/clair.yml'
	shell:
		'vcfcat {input.ct1} {input.ct2} > {output}'
rule prep4Bcf:
	input:
		rules.vcfcat.output
	output:
		'{results}/{sample}_flye_clair.norm.vcf.gz'
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/clair.yml'
	shell:
		'bgzip {input} && tabix -p vcf {output}'

rule bcftools:
	input:
		assembly = rules.flye.output,
		bcfprep = rules.prep4Bcf.output
	output:
		'{results}/{sample}_flye.clair.fasta'
	conda:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/clair.yml'
	shell:
		'bcftools consensus -f {input.assembly} {input.bcfprep} > {output}'
rule clairPolca:
	input:
		gen = rules.bcftools.output,
		r1 = expand('{sampleDir}/{r1}',sampleDir=sampleDir,r1=config["illumina"][0]),
                r2 = expand('{sampleDir}/{r2}',sampleDir=sampleDir,r2=config["illumina"][1])
	output:
		'{results}/{sample}ClairPolca'
	shell:
		"polca.sh -a {input.gen} -r '{input.r1} {input.r2}' && mkdir {output} && mv {wildcards.sample}_flye.clair.fasta* {output}"	
	
