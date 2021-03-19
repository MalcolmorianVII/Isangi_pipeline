#report:"report/workflow.rst"	# Global description of the workflow
configfile:"config.yaml"
#samples = ['CHF10J']
root_dir = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.02.06']
res=['/home/ubuntu/data/belson/isangi_nanopore/qc/results/pipeline_optimization/flye_polca_medaka']
rule all:
	input:
		expand('{res}/{sample}_polish_medaka_results',res=res,sample=config["nanopore"]),
		expand('{res}/{sample}_polish_flye_results',res=res,sample=config["nanopore"])

rule denovo:
	input:
		expand('{root}/{sample}.fastq.gz',root=root_dir,sample=config["nanopore"])
	output:
		report(directory('{res}/{sample}_flye'), patterns=["assembly.fasta"],caption="report/denovo.rst", category="Step 1")
	shell:
		'flye --nano-raw {input} -g 5m -o {output} -t 8 --plasmids'

rule polish_flye:
	input:
		assembly = rules.denovo.output,
		r1 = expand('{root}/{r1}',root=root_dir,r1=config["illumina"][0]),
		r2 = expand('{root}/{r2}',root=root_dir,r2=config["illumina"][1])
	output:
		report(directory('{res}/{sample}_polish_flye_results'),patterns=["assembly.fasta.PolcaCorrected.fa"],caption="report/polish_flye.rst",category="Step 2")
	shell:
		"polca.sh -a {input.assembly}/assembly.fasta -r '{input.r1} {input.r2}' && mkdir {output} && mv assembly.fasta* {output}"

rule medaka:
        input:
                rules.polish_flye.output
        output:
               report(directory('{res}/{sample}_medaka'),patterns=["consensus.fasta"],caption="report/medaka.rst",category="Step 3")
        conda:
                'envs/medaka.yml'
        shell:
                'medaka_consensus -i {rules.denovo.input} -d {input}/assembly.fasta.PolcaCorrected.fa -t 8  -m r941_min_high_g360 -o {output}'

rule polish_medaka:
        input:
                rules.medaka.output
        output:
                report(directory('{res}/{sample}_polish_medaka_results'),patterns=["consensus.fasta.PolcaCorrected.fa"],caption="report/polish_medaka.rst",category="Step 4")
        shell:
                "polca.sh -a {input}/consensus.fasta -r '{rules.polish_flye.input.r1} {rules.polish_flye.input.r2}' && mkdir {output} && mv consensus.fasta* {output}"
