#! /usr/bin/bash

PolcaAllResults=${1}Results.txt
for dir in ${1}/*olish*
do
	echo $dir >> ${1}/$PolcaAllResults
	cat $dir/*report >> ${1}/$PolcaAllResults
	
done
echo "Malcolm 7 Results aggrgation done"
