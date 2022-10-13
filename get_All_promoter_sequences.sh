
#!/bin/bash

#Script used to excract promoter regions (1kb) of all genes within the GCN
#Necessary inputs/files:
# - annotation file with all cork oak genes from CorkOakDB webiste: GCF_002906115.1_CorkOak1.0_genomic.gff > CorkOak_genes.gff
#Run command used: getPromoterRegions.sh ./GCF_002906115.1_CorkOak1.0_genomic.gff 1000

#Getting a bed file with all cork oak genes (without trnas and mrnas)
grep -P '\tgene\t' $1 > ./CorkOak_genes.gff
awk '{print $1,$4,$5,$9,$7}' ./CorkOak_genes.gff > ./New_CorkOak_genes.gff
sed 's/ID=gene-//g' ./New_CorkOak_genes.gff > ./New_CorkOak_genes2.gff
sed 's/;.* / /g' ./New_CorkOak_genes2.gff > ./New_CorkOak_genes3.bed
sed 's/ /\t/g' ./New_CorkOak_genes3.bed > ./CorkOak_genes.bed

#Getting chromosome/scafold sizes in a fasta file
grep -P '\tregion\t' $1 > ./CorkOak_scafolds.gff
awk '{print $1,$5}' ./CorkOak_scafolds.gff > ./CorkOak_scafolds.fa
sed 's/ /\t/g' ./CorkOak_scafolds.fa > ./CorkOak_scafold_sizes.fa

#Getting promoter regions of all genes with "$2" kb
grep '+' ./CorkOak_genes.bed > ./CorkOak_genes_p.bed
bedtoolds flank -i ./CorkOak_genes_p.bed -g ./CorkOak_scafold_sizes.fa -l $2 -r 0 > ./CorkOak_promoters_1kb.bed
grep '-' ./CorkOak_genes.bed > ./CorkOak_genes_n.bed
bedtools flank -i ./CorkOak_genes_n.bed -g ./CorkOak_scafold_sizes.fa -l 0 -r $2 >> ./CorkOak_promoters_1kb.bed

sed 's/ .*//g' $1 > ./CorkOak_scafolds.fna
bedtools getfasta -fi ./CorkOak_scafolds.fna -bed ./CorkOak_promoters_1kb.bed -fo ./CorkOak_promoter_1kb_sequences.bed

