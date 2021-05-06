#! /usr/bin/bash
set -e

D=`date +"%Y.%m.%d"`
A=/home/ubuntu/data/belson/isangi_nanopore/qc/results/analysis
P=/home/ubuntu/data/belson/isangi_nanopore/qc/results/polishing
 
echo "Choose destination directory:a-Analysis & p-Polishing:"

read REPLY

if [ $REPLY == "a" ]
then
	mkdir -p $A/$D | ls $A

elif  [ $REPLY == "p" ]
then
	mkdir -p $P/$D | ls $P

else
	echo "Invalid argument"
fi

exit 0
