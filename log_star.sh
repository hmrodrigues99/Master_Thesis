#!/bin/bash

#This script receives the SRR_pass_filtered.Log.final.out files from the STAR run and creates a STAR_log.csv with the number of uniquely mapped reads and its percentage in the second and third collumns. The fourth collumn shows the sum of multi-mapping, unmapped and chimeric reads. The fifth collumn only shows the number of unmapped reads (due to lack of annotation knowledge). The first collumn corresponds to the sample name, one per line.

echo $1

#number of input reads
input_reads_n=$(grep "Number of input reads" $1 | grep -o '[[:digit:]]*')
echo $input_reads_n

#uniquely mapped reads number
uniq_map_reads_n=$(grep "Uniquely mapped reads number" $1 | grep -o '[[:digit:]]*')
echo $uniq_map_reads_n

#uniquely mapped reads percentage
uniq_map_reads_p=$(grep "Uniquely mapped reads %" $1 | grep -o '[[:digit:]][[:digit:]]*.[[:digit:]]*')
echo $uniq_map_reads_p

#total multi-mapped/unmapped and chimeric reads
mm_um_cr_n=$(($input_reads_n-$uniq_map_reads_n))
echo $mm_um_cr_n

#number of unmapped reads (due to lack of annotation)
unmapped_reads_short=$(grep "Number of reads unmapped: too short" $1 | grep -o '[[:digit:]]*')
unmapped_reads_other=$(grep "Number of reads unmapped: other" $1 | grep -o '[[:digit:]]*')
unmapped_reads_n=$(($unmapped_reads_short+$unmapped_reads_other))
echo $unmapped_reads_n

printf "${uniq_map_reads_n} \t ${uniq_map_reads_p} \t ${mm_um_cr_n} \t ${unmapped_reads_n} \n" >> STAR_log.csv
