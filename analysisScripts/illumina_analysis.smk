configfile:"config.yaml"
#todo_list = ['17761-33892_1_65','17762-33892_1_71','17763-33892_1_73','17764-33892_1_77']
root_dir = '/home/ubuntu/samples'

#rule_order: Denovo > mlst,assemblt-stats

rule all:
	input:
		mlst = expand('/home/ubuntu/samples/{sample}/{sample}_mlst.txt',sample = config["samples"]),
		assem = expand('/home/ubuntu/samples/{sample}/{sample}_assembly.txt',sample = config["samples"]),
		report = expand('{sample}_report.txt',sample = config["samples"])
rule denovo:
	input:
		r1 = '/home/ubuntu/samples/{sample}_1.fastq',
		r2 = '/home/ubuntu/samples/{sample}_2.fastq'
	output:
		out = directory('/home/ubuntu/samples/shovill/{sample}')
	shell:
		'shovill --outdir /home/ubuntu/samples/shovill/{wildcards.sample} --R1 {input.r1} --R2 {input.r2} --cpus 4 --ram 16'


rule mlst:
        input:
                rules.denovo.output.out
        output:
                '/home/ubuntu/samples/{sample}/{sample}_mlst.txt'
        shell:
                'mlst --scheme senterica --nopath {input}/contigs.fa > {output}'

rule assembly:
        input:
                rules.denovo.output.out
        output:
                '/home/ubuntu/samples/{sample}/{sample}_assembly.txt'
        shell:
                'assembly-stats {input}/contigs.fasta > {output}'
rule done:
	input:
		rules.mlst.output,
		rules.assembly.output
	output:
		'{sample}_report.txt'
	script:
		'scripts/report.py'
