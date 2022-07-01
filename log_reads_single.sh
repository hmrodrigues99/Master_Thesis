#!/bin/bash

#This script receives a trimmomatic.log and outputs a trimmomatic_quality.csv  with the number of kept and discarded numbers of reads, and the survival rate of reads for each sample (one line for sample).

#Total number of reads
total_reads_prov=$(wc -l $1)
total_reads=$(echo $total_reads_prov | sed 's/ .*//')
echo $total_reads

#Number of discarded single 0 0 0 reads
discarded_reads=$(grep ".*[[:space:]]0[[:space:]]0[[:space:]]0.*" $1 | wc -l)
echo $discarded_reads

#Survival rate of reads
survival_rate_prov=$(($total_reads - $discarded_reads))
survival_rate_prov2=$(($survival_rate_prov * 100))
survival_rate=$(echo "scale=5 ; $survival_rate_prov2 / $total_reads" | bc)
echo $survival_rate

printf "$1 \t $total_reads \t $discarded_reads \t $survival_rate \n" >> trimmomatic_log.csv

