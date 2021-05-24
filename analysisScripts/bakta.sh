#!/usr/bin/bash

#NOTE:  $1 = {input}; $2 =  {output}
set -e
eval "$(conda shell.bash hook)"
conda activate "bakta"

bakta --db /home/ubuntu/data/belson/bakta/db/db $1 -o $2
exit 0
