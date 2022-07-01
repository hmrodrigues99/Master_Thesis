
#!/bin/bash

#This script receives the BioProject_raw_counts.txt.summary and outputs a featureCounts_log.csv where each sample corresponds to a line, first collumn the Assigned reads and the second collumn lists the Unassigned_NoFeatures

printf "Assigned Reads \n" > feat1.csv
printf "Unassigned Reads \n" > feat2.csv

#Number of Assigned Reads
ass_reads_n=$(grep "Assigned" $1 | grep -o '[[:digit:]]*')

#Number of Unassigned_NoFeatures
unass_reads_n=$(grep "Unassigned_NoFeatures" $1 | grep -o '[[:digit:]]*')

for i in $ass_reads_n;do printf "$i \n" >> feat1.csv; done
for j in $unass_reads_n;do printf "$j \n" >> feat2.csv; done
paste feat1.csv feat2.csv > featureCounts_log.csv

rm $PWD/feat1.csv
rm $PWD/feat2.csv
