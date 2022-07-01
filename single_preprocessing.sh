
#!/bin/bash

#This script receives a list of samples $1 SINGLE_END READS and runs FastQC, Trimmomatic 0-39 and FastQC again, creating a preprocessing folder on the BioProject directory in which the script was ran.

export PATH=$PATH:~/bin/FastQC # in server might be different
export PATH=$PATH:~/data/scripts # path to the log_reads script
trim_path="$HOME/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters="$HOME/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
mkdir -p preprocessing/fastqc/raw
mkdir -p preprocessing/fastqc/filtered
mkdir -p preprocessing/hq_reads
printf "Sample \t Total Reads \t Discarded Reads \t Trimmed Reads (perc) \n" > trimmomatic_log.csv

#First run of FastQC for the raw reads

list=""
for f in $(ls -1 fastq/*fastq.gz)
do
        list+=" ${f}"
done

echo "Starting FastQC for raw reads"
#fastqc -t 6 -o preprocessing/fastqc/raw $list

#Running trimmomatic-0.39.jar
for f in $(ls -1 fastq/*fastq.gz)
do
        sample=$(basename $f | sed "s/.fastq.gz//g") 
        echo "Starting Trimmomatic for ${sample}"
        echo "java -jar $trim_path SE $f -baseout preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
        java -jar $trim_path SE $f preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        #Running log_reads_single.sh for each sample.log
        log_reads_single.sh preprocessing/hq_reads/$sample.filtered.log $sample

done

#Second run of FastQC for the filtered reads

list=""
for f in $(ls -1 preprocessing/hq_reads/*.fq.gz)
do
        list+=" ${f}"
done

echo "Starting FastQC for filtered read"
fastqc -t 6 -o preprocessing/fastqc/filtered $list
