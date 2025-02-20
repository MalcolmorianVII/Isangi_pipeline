rule all:
    input:
        "/data/fast/core/isangi/results/CHF10J1.guppy5.sup_racon.fasta",
        "/data/fast/core/isangi/results/CHF10J1.guppy5.sup_medaka",
        # "/data/fast/core/isangi/results/CHF10J1.guppy5.sup_circlator",
        "/data/fast/core/isangi/results/amr_finder_plus/CHF10J1.guppy5.sup.amr_finder_plus.tsv",
        "/data/fast/core/isangi/results/CHF10J1.guppy5.sup_bakta"


rule denovo:
    # First assembling the genome using Flye
    input:
        "/data/fast/core/isangi/{sample}.fastq.gz",
    output:
        directory("/data/fast/core/isangi/results/{sample}_flye"),
    conda:
        "/home/bkutambe/miniconda3/envs/flye"
    shell:
        "flye --nano-hq {input} -g 5m -o {output} -t 32"


# Map the accurate short reads to the assembly for polishing
rule minimap:
    input:
        ass="/data/fast/core/isangi/results/{sample}_flye/assembly.fasta",
        r1="/data/fast/core/isangi/17762-33892_1_71_bbduk_1.fastq.gz",
        r2="/data/fast/core/isangi/17762-33892_1_71_bbduk_2.fastq.gz",
    output:
        "/data/fast/core/isangi/results/{sample}.paf",
    conda:
        "/home/bkutambe/miniconda3/envs/minimap"
    shell:
        "minimap2 -x sr {input.ass} {input.r1} {input.r2} > {output}"


# Genome polishing with Racon
rule racon:
    input:
        paf=rules.minimap.output,
        r1="/data/fast/core/isangi/17762-33892_1_71_bbduk_1.fastq.gz",
        #r2 = "/data/fast/core/isangi/17762-33892_1_71_bbduk_2.fastq.gz",
        ass="/data/fast/core/isangi/results/{sample}_flye/assembly.fasta",
    output:
        "/data/fast/core/isangi/results/{sample}_racon.fasta",
    conda:
        "/home/bkutambe/miniconda3/envs/racon"
    shell:
        "racon -t 4 {input.r1} {input.paf} {input.ass} > {output}"


#  Genome polishing with Medaka
rule medaka:
    input:
        racon=rules.racon.output,
        read=rules.denovo.input
    output:
        directory("/data/fast/core/isangi/results/{sample}_medaka")
    conda:
        "/home/bkutambe/miniconda3/envs/medaka"
    shell:
        "medaka_consensus -i {input.read} -d {input.racon} -t 8  -m r941_min_sup_g507 -o {output}"


# Circularizing the draft genome after polishing i.e correcting breaks introduced during polishing & assembly????
# rule circlator:
#     input:
#         medaka = rules.medaka.output,
#         read=rules.denovo.input
#     output:
#         directory("/data/fast/core/isangi/results/{sample}_circlator")
#     container:
#         "docker://sangerpathogens/circlator:latest"
#     shell:
#         "circlator all --merge_min_id 85 --merge_breaklen 1000 {input.medaka}/consensus.fasta {input.read} {output}"

rule amr_finder_plus:
   input:
       assembly = rules.medaka.output
   output:
       amr_finder_plus_results = '/data/fast/core/isangi/results/amr_finder_plus/{sample}.amr_finder_plus.tsv'
   conda:
       '/home/bkutambe/miniconda3/envs/amrfinder'
   shell:
       'amrfinder -n {input.assembly}/consensus.fasta -O Salmonella --output {output.amr_finder_plus_results} --threads 4 --name {wildcards.sample}'

rule bakta:
   input:
       assembly = rules.medaka.output
   output:
       directory("/data/fast/core/isangi/results/{sample}_bakta")
   conda:
       '/home/bkutambe/miniconda3/envs/bakta'
   shell:
       'bakta --db /data/fast/core/bakta/db --threads 16 --output {output} --prefix {wildcards.sample} --force {input.assembly}/consensus.fasta'
