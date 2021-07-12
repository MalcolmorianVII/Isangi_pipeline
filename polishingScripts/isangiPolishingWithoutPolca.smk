# Genomic analysis workflow i.e denovo assembly + genome polishing + annotation
#configfile:"/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/isangiConfig.yaml"
barcode = ["barcode01","barcode02","barcode03","barcode04","barcode05","barcode06","barcode07","barcode08","barcode09","barcode10","barcode11","barcode12"]
root = "/home/ubuntu/data/raw_data/20201203_1331_MN33881_FAO37426_b676e9be_guppy5_sup/fastq_pass"
#root_dir = ['/home/ubuntu/data/belson/test_data/2021.01.04']
long = "fastq_runid_96ecb7a94086f36bb155ec67c59e89d445e1275d_0.fastq.gz"
results = ['/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing/2021.07.08']
#nano = '/home/ubuntu/data/belson/monash_data/guppy_v3.6.1_reads'
def nano(wildcards):
	 return expand('{root}/{sample}/{long}',root=root,sample=wildcards.sample,long=long)
rule all:
	input:
		expand('{results}/{sample}/{sample}Flye',sample = barcode,results = results),
		expand('{results}/{sample}/{sample}RaconX1.fasta',sample = barcode,results = results),
		#expand('{results}/{sample}RaconX2.fasta',sample = barcode,results = results),
		#expand('{results}/{sample}RaconX3.fasta',sample = barcode,results = results),
		#expand('{results}/{sample}RaconX4.fasta',sample = barcode,results = results),
		expand('{results}/{sample}/{sample}Medaka',sample = barcode,results = results)

rule flye:
	input:
		L = nano
		#expand('{root}/{{sample}}/{nano}',root=root,sample=barcode,nano=nano)
	output:
                directory('{results}/{sample}/{sample}Flye')
	shell:
                'flye --nano-raw {input.L} -g 5m -o {output} -t 8 --plasmids'

rule raconX1:
        input:
                gen = rules.flye.output,
		L = nano
        output:
                x1 = '{results}/{sample}/{sample}RaconX1.fasta',
                pf1 = temp('{results}/{sample}.racon.paf')
        shell:
                'minimap2 -x map-ont {input.gen}/assembly.fasta {input.L} > {output.pf1} && racon -t 4 {input.L} {output.pf1} {input.gen}/assembly.fasta > {output.x1}'


rule medaka:
        input:
                gen = rules.raconX1.output.x1,
		L = nano
        output:
                directory('{results}/{sample}/{sample}Medaka')
        conda:
                '/home/ubuntu/data/belson/isangi_nanopore/qc/scripts/envs/medaka.yml'
        shell:
                'medaka_consensus -i {input.L} -d {input.gen} -t 8  -m r941_min_sup_g507 -o {output}'
