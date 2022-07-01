#!/bin/bash

#This script receives as input $1 a Bioprojects.txt with each Illumina Bioproject name and runs for each Bioproject the multi_fastqdump.sh which retrieves the SRA sequences for each one of them.

while read f
do
	$HOME/data/scripts/multi_fastqdump.sh $f
done < $1
