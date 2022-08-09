
#!/bin/bash

#This script runs the STAR software with the filtered sequences located in hq_reads. It begins with a unique run to generate a genome_index after the genome indexes are generated
#Should be run in the specific BioProject directory
#FOR SINGLE-END READS ONLY!

export PATH=$PATH:$HOME/bin/STAR-2.7.9a/bin/Linux_x86_64
export PATH=$PATH:~/data/scripts
mkdir -p aligned_reads
gindex="$HOME/data/genome/genome_index"
grindex="$HOME/data/genome"
bioproject_name=$(echo $PWD | sed "s,.*data/,,g")
printf "Unique Mapped Reads \t Uniq Mapped Reads (perc) \t Total Unmapped Reads \t Unmapped Reads (Anotation) \n" > STAR_log.csv

#Running STAR to create genome indexes

#STAR --genomeDir genomeGenerate --runThreadN 10 --genomeDir $gindex --genomeFastaFiles $grindex --sjdbGTFfile "path to annotations.gff" --sjdbGTFtagExonParentTranscript Parent

#Running STAR with the samples in hq_reads
for f in $(ls -1 $PWD/preprocessing/hq_reads/*.fq.gz)
do
        sample=$(basename $f)
        echo "Starting STAR for $sample"
        STAR --runThreadN 20 --genomeDir $gindex --readFilesIn $f --outFileNamePrefix aligned_reads/$sample. --twopassMode Basic --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 100000 --outSAMtype BAM Unsorted

	#Gathering the data from the Log.final.out files from all samples in the Bioproject and making a STAR_log.csv
	log_star.sh aligned_reads/$sample.Log.final.out
done

#This script receives a list of BAM files and outputs a bioproject_nameraw_counts.txt.
#Should be run in the specific Bioproject directory.

export PATH=$PATH:~/bin/subread-2.0.3-Linux-x86_64/bin

bam_list=""
for i in $(ls -1 aligned_reads/*Aligned.out.bam); do bam_list+=" ${i}"; done

featureCounts -T 5 -t exon -g gene -a ~/data/genome/GCF_002906115.1_CorkOak1.0_genomic.gff -o ${bioproject_name}raw_counts.txt $bam_list

#This script uses the raw_counts.txt.summary to compute a provisional featureCounts_log.csv
log_featureCounts.sh ${bioproject_name}raw_counts.txt.summary 's'
#This script merges the trimmomatic_quality.csv, the STAR_log.csv and the featureCounts_log.csv into a BioProject_Statistics.csv
log_final.sh ${bioproject_name}
