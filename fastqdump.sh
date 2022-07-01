#!/bin/bash

export PATH=$PATH:~/bin/sratoolkit.2.11.1-ubuntu64/bin/ #path to fastq-dump

#Input $1 bioproject.txt that contains the SRA sequences.
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
