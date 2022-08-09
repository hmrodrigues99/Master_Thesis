
#!/bin/bash

#This script receives the $1 BioProject_raw_counts.txt.summary and a $2 signaling string: 's' for single reads and 'p' for paired end reads.
It outputs a featureCounts_log.csv where each line represents a sample, first collumn counts Assigned reads and the second collumn counts the Unassigned_NoFeatures reads.

read_type=$2
printf "Assigned Reads \n" > feat1.csv
printf "Unassigned Reads \n" > feat2.csv

#Number of Assigned Reads
ass_reads_n=$(grep "Assigned" $1 | grep -o '[[:digit:]]*')

#Number of Unassigned_NoFeatures
unass_reads_n=$(grep "Unassigned_NoFeatures" $1 | grep -o '[[:digit:]]*')

#Check if the reads are paired, and if they are, double the read counts
if [ $read_type == 'p' ];then
  ass_reads_double_n2=$((ass_reads_n))
  ass_reads_double_n=$(($ass_reads_double_n2*2))
  unass_reads_double_n2=$((unass_reads_n))
  unass_reads_double_n=$((unass_reads_double_n2*2))
  
  for i in $ass_reads_double_n;do printf "$i \n" >> feat1.csv; done
  for j in $unass_reads_double_n;do printf "$j \n" >> feat2.csv; done
else
  for i in $ass_reads_n;do printf "$i \n" >> feat1.csv; done
  for j in $unass_reads_n;do printf "$j \n" >> feat2.csv; done
fi

paste feat1.csv feat2.csv > featureCounts_log.csv

rm $PWD/feat1.csv
rm $PWD/feat2.csv
