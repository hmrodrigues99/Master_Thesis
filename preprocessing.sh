#!/bin/bash

#This script receives a list of samples $1, the sufix for the first pair $2 and the sufix for the second pair in $3

export PATH=$PATH:~/bin/FastQC # in server might be different
export PATH=$PATH:~/data/scripts # path to the log_reads script
trim_path="$HOME/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters="$HOME/bin/Trimmomatic-0.39/adapters/TruSeq2-PE.fa"
mkdir -p preprocessing/fastqc/raw
mkdir -p preprocessing/fastqc/filtered
mkdir -p preprocessing/hq_reads
r1_suf=$2
r2_suf=$3
printf "Sample \t Total Reads \t Discarded Reads \t Trimmed Reads (perc) \n" > trimmomatic_log.csv

#First run of FastQC for the raw reads

while read f
do
	in1="${f}${r1_suf}"
	in2="${f}${r2_suf}"
        date
	sample=$(basename $f)
	echo "Starting FastQC for ${sample}"
	echo "fastqc -t 6 -o preprocessing/fastqc/raw ${in1} ${in2}"
        fastqc -t 6 -o preprocessing/fastqc/raw ${in1} ${in2}
		

	#Running trimmomatic-0.39.jar	
	echo "Starting Trimmomatic for ${sample}"
	echo "java -jar $trim_path PE $in1 $in2 -baseout preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
	java -jar $trim_path PE $in1 $in2 -baseout preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	#Running log_reads.sh for each sample.log
	log_reads.sh preprocessing/hq_reads/$sample.filtered.log

done < $1

#Second run of FastQC for the filtered reads

list=""
for f in $(ls -1 preprocessing/hq_reads/*P.fq.gz)
do
	list+=" ${f}"
done    

echo "Starting FastQC for filtered read"
fastqc -t 6 -o preprocessing/fastqc/filtered $list
