#! /usr/bin/bash
#set -e
set +eu
eval "$(conda shell.bash hook)"
conda activate snakemake
# workflowgraph dir
wf=/home/ubuntu/data/belson/bioinformatics/projects_2021/isangi_nanopore/scripts/graphs

# $1=script

base=`basename -s .smk $1`
wfg=$wf/$base.svg
snakemake --dag --snakefile $1 | dot -Tsvg > $wfg 
if [ -f $wfg ]
then
	echo $wfg created!!!!
else
	echo an error occured!
fi
exit 0
