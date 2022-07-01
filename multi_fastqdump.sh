#!/bin/bash

'''
General Use
This script, in conjunction with multi_fastqdump.sh, performs the download from NCBI of RNA-seq libraries by their BioProject accession ID and SRA run codes.
Inputs:
1) A list of BioProject accession IDs stored into NCBI (ex: PRJNA392919), one ID per line.
Outputs:
For each BioProject, a new directory is created to store the retrieved data (ex: PRJNA392919/fastq).
For each BioProject, a file with the retrieved read names is also created (ex: PRJNA392919/PRJNA392919_raw_sample_list).

It is necessary to have in the same run directory of this script, a number of lists, one for each BioProject, where the name of the list is
the BioProject accession ID and the contents are the SRA run codes for that BioProject, one per line (ex: SRR5986737). 
The SRA run codes need to be manually picked for each correspondent NCBI BioProject.

Thesis Use
fastq_dump.sh Bioprojects.txt
This script received as input $1 a file named Bioprojects.txt with each Illumina Bioproject accession ID.
For each BioProject there was a list named with the BioProject accession ID, and contained the SRA run codes for that BioProject. (ex:PRJNA392919.txt)
'''

while read f
do
	#Script running directory
	$HOME/data/scripts/multi_fastqdump.sh $f
done < $1
