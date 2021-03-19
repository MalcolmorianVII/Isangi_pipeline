configfile:"config.yaml"
#samples = ['CHF10J']
root_dir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
res=['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.03.08']
rule all:
	input:
		expand('{res}/{sample}_plasmid_amrfinder',res=res,sample=config["nanopore"]),
		expand('{res}/{sample}_plasmid_prokka',res=res,sample=config["nanopore"]),
		expand('{res}/{sample}_chromosome_amrfinder',res=res,sample=config["nanopore"]),
		expand('{res}/{sample}_chromosome_prokka',res=res,sample=config["nanopore"])

rule denovo:
	input:
		expand('{root}/{sample}.fastq.gz',root=root_dir,sample=config["nanopore"])
	output:
		directory('{res}/{sample}_flye')
	shell:
		'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polish_flye:
	input:
		assembly = rules.denovo.output,
		r1 = expand('{root}/{r1}',root=root_dir,r1=config["illumina"][0]),
		r2 = expand('{root}/{r2}',root=root_dir,r2=config["illumina"][1])
	output:
		directory('{res}/{sample}_polish_flye_results')
	shell:
		"polca.sh -a {input.assembly}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule plasmid:
	input:
		rules.polish_flye.output
	output:
		'{res}/{sample}_isangi_plasmid.fasta'
	shell:
		"./plasmid_genome_splitter.sh plasmid {input}/assembly.fasta.PolcaCorrected.fa {output}"
rule chromosome:
	input:
		rules.polish_flye.output
	output:
		'{res}/{sample}_chromosome.fasta'
	shell:
		"./plasmid_genome_splitter.sh chromosome {input}/assembly.fasta.PolcaCorrected.fa {output}"		
rule amrfinder:
	input:
		plasmid=rules.plasmid.output,
		chromo = rules.chromosome.output

	output:
		plasmid = '{res}/{sample}_plasmid_amrfinder',
		chromo = '{res}/{sample}_chromosome_amrfinder'
	shell:
		"""
		amrfinder --plus -n {input.plasmid} -O Salmonella > {output.plasmid}
		amrfinder --plus -n {input.chromo} -O Salmonella > {output.chromo}
		"""

rule prokka:
	input:
		plasmid = rules.plasmid.output,
		chromo = rules.chromosome.output
	output:
		plasmid = directory('{res}/{sample}_plasmid_prokka'),
		chromo = directory('{res}/{sample}_chromosome_prokka')
	conda:
		'envs/prokka.yml'
	shell:
		"""
		prokka {input.plasmid} --outdir {output.plasmid}
		prokka {input.plasmid} --outdir {output.chromo}
		"""
