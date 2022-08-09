
#!/bin/bash

#This script runs the STAR software with the filtered sequences located in hq_reads. It begins with a unique run to generate a genome_index after the genome indexes are generated
#Should be run in the specific BioProject directory
#FOR PAIRED-END READS ONLY!

export PATH=$PATH:$HOME/bin/STAR-2.7.9a/bin/Linux_x86_64
export PATH=$PATH:~/data/scripts
mkdir -p aligned_reads
gindex="$HOME/data/genome/genome_index"
grindex="$HOME/data/genome"
bioproject_name=$(echo $PWD | sed "s,.*data/,,g")
printf "Unique Mapped Reads \t Uniq Mapped Reads (perc) \t Total Unmapped Reads \t Unmapped Reads (Anotation) \n" > STAR_log.csv

#Running STAR to create genome indexes

#STAR --genomeDir genomeGenerate --runThreadN 10 --genomeDir $gindex --genomeFastaFiles $grindex --sjdbGTFfile "path to annotations.gff" --sjdbGTFtagExonParentTranscript Parent

#First it creates a samples.txt with the readsP from all samples
ls -1 $PWD/preprocessing/hq_reads/*P.fq.gz | sed "s/_.P.*//g" | uniq | sort > $PWD/filtered_samples.txt
r1_suf="_1P.fq.gz"
r2_suf="_2P.fq.gz"
#Running STAR with the samples
while read -r line
do
	in1="${line}${r1_suf}"
	in2="${line}${r2_suf}"
	sample=$(basename $line)
	echo "Starting STAR for $sample"
	STAR --runThreadN 20 --genomeDir $gindex --readFilesIn $in1 $in2 --outFileNamePrefix aligned_reads/$sample. --twopassMode Basic --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 100000 --outSAMtype BAM Unsorted

	#Gathering the data from the Log.final.out files from all samples in the Bioproject and making a STAR_log.csv
	log_star.sh aligned_reads/$sample.Log.final.out

done < $PWD/filtered_samples.txt

#The following code receives a list of BAM files and outputs a bioproject_nameraw_counts.txt.
export PATH=$PATH:$HOME/bin/subread-2.0.3-Linux-x86_64/bin

#Sort the BAM file with samtools sort
samtools sort -o "output sorted bam" -@ 5 "input unsorted bam"

bam_list=""
for i in $(ls -1 aligned_reads/*Aligned.out.bam); do bam_list+=" ${i}"; done

featureCounts -T 5 -p --countReadPairs -t exon -g gene -a ~/data/genome/GCF_002906115.1_CorkOak1.0_genomic.gff -o ${bioproject_name}raw_counts.txt $bam_list

#This script uses the raw_counts.txt.summary to compute a provisional featureCounts_log.csv
log_featureCounts.sh ${bioproject_name}raw_counts.txt.summary 'p'
#This script merges the trimmomatic_quality.csv, the STAR_log.csv and the featureCounts_log.csv into a BioProject_Statistics.csv
log_final.sh ${bioproject_name}
