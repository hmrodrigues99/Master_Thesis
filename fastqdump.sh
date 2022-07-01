#!/bin/bash

'''
This script runs in loop under the multi_fastqdump.sh.
This script can be ran alone to retrieve the RNA-seq libraries of a single BioProject accession ID through the SRA (sequence read archive) from NCBI.
Input:
1)File named with the BioProject accession ID (ex. PRJNA392919.txt), containing the SRA run codes correspondent to that BioProject, one per line (ex: SRR5986737).

#Path to the program fastq-dump
export PATH=$PATH:~/bin/sratoolkit.2.11.1-ubuntu64/bin/

#Input $1 Bioproject.txt that contains the SRA sequences.
#Makes new directory, if it doesnt exist already, with the Bioproject access ID.
dirname=$(echo "$1" | cut -f 1 -d '.')
if [ ! -d $dirname ]; then
	mkdir -p $dirname
fi

while read f
do
	echo $f
	date
	echo "fastq-dump --outdir $dirname/fastq --gzip --skip-technical --readids --read-Filter pass --dumpbase --split-3 --clip $f"
	fastq-dump --outdir $dirname/fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $f
	echo "DONE"
done < $1

ls -1 ~/data/$dirname/fastq/*_1.fastq.gz | sed "s/_1.fastq.gz//g" > ~/data/$dirname/${dirname}_raw_sample_list.txt
