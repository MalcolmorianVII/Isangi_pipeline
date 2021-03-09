#!/usr/bin/bash
# $1 specifies what output to generate i.e plasmid or chromosome
# $2 the input file i.e assembly.fasta in this case
# $3 the output i.e plasmid.fasta or chromosome.fasta according to $1
set +e
num=$(grep -n "contig_1" $2 | cut -d: -f1)

if  [ $1 == "plasmid" ]	#For plasmids
then
	num=$(expr $num - 1)
	head -n $num $2 > $3
	echo Created $3 file

elif [ $1 == "chromosome" ] 	#For chromosomes
then
	num=$( expr $(wc -l assembly.fasta.PolcaCorrected.fa | cut -d " " -f1) - $num + 1)
	tail -n $num $2 > $3
	echo Created $3 file
	#echo $num

else
	echo Incorrect output specified!!!
fi
echo Done!!!!!!
