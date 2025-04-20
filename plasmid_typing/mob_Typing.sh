#! /usr/bin/bash
#url=http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&sendto=on&id=LC542971.1

while read line 
do
	mkdir -p /home/ubuntu/data/belson/mob-suite/mob_results/$line
	wget -O /home/ubuntu/data/belson/mob-suite/mob_results/$line/$line.gb "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&dopt=fasta&sendto=on&id=$line"
	mob_typer --infile=/home/ubuntu/data/belson/mob-suite/mob_results/$line/$line.gb --out_file=/home/ubuntu/data/belson/mob-suite/mob_results/$line/${line}Mobtyper.txt 
done < $1
echo "MALCOLM program finished with 0 exit code"
