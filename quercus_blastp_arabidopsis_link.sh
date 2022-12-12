
#!/bin/bash

#This script performs a blastp to obtain a sequence homology link between Quercus suber and Arabidopsis thaliana genes

#Necessary commands installed:
# - makeblastdb
# - blastp
#Necessaru inputs/files:
# - arabidopsis proteome fasta file from TAIR website: Araport11_pep_20220193.faa
# - quercus suber proteome fasta file from CorkOakDB website: suber_proteome
# - name of the output file: quercus_blastp_thaliana_top20.txt

#Run command used: quercus_blastp_arabidopsis.sh Araport11_pep_20220193.faa suber_proteome.fasta quercus_blastp_thaliana_top20.txt

#Creating the arabidopsis database
makeblastdb - query $1 -db arabidopsis_thaliana_DB - dbtype prot -parse_seqids

#Running the blastp
blastp -query $2 -db arabidopsis_thaliana_DB -max_target_seqs 20 -evalue 1e6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue' -out quercus_blastp_thaliana_top20.txt

