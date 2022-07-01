#!/bin/bash

#This script receives a trimmomatic.log and outputs a trimmomatic_quality.csv  with the number of kept and discarded numbers of reads, and the survival rate of reads for each sample (one line for sample).

#Total number of reads
total_reads_prov=$(wc -l $1)
total_reads=$(echo $total_reads_prov | sed 's/ .*//')

#Number of discarded single 0 0 0 reads
discarded_single_reads=$(grep ".*[[:space:]]0[[:space:]]0[[:space:]]0.*" $1 | wc -l)

#Number of pairs constituted by both 0 0 0 reads
pairs=$(pcregrep -Mc '.*1[[:space:]]0[[:space:]]0[[:space:]]0[[:space:]]0\n.*2[[:space:]]0[[:space:]]0[[:space:]]0' $1)

#Number of discarded reads
pairs2=$((pairs * 2))
other_pairs=$(($discarded_single_reads - $pairs2))
other_pairs2=$((other_pairs * 2))
discarded_reads=$((pairs2 + other_pairs2))

#Survival Rate
survival_rate_prov=$(($total_reads - $discarded_reads))
survival_rate_prov2=$(($survival_rate_prov * 100))
survival_rate=$(echo "scale=5 ; $survival_rate_prov2 / $total_reads" | bc)

#Sample Name
sample=$(echo $1 | sed "s,.*reads/,,g" | sed "s,_pass.*,,g")

printf "$sample \t $total_reads \t $discarded_reads \t $survival_rate \n" >> trimmomatic_log.csv 
