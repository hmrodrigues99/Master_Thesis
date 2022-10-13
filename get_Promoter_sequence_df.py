# -*- coding: utf-8 -*-

#This script creates a table (Promoter_sequences_df.csv) with the cork oak genes and their respective promoter sequences that are present in the Cork network.

#Necessary inputs/files
# - Cork network: Cork_network.txt
# - Promoter regions of all cork oak genes: CorkOak_promoter_1kb_sequences.bed)
# - Promoter region coordinates of all cork oak genes: CorkOak_promoters_1kb.bed

import pandas as pd

#input Cork network
network = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_2/Cork_network.txt")

#Input file with all the gene promoter sequences
promoter_sequences = pd.read_table("C:/Users/Hugo Rodrigues/Documents/TF_Targets/Promoters/CorkOak_promoter_1kb_sequences.bed", sep="\t", names=["Scafolds"])

#Input promoter coordinates
promoter_co = pd.read_table("C:/Users/Hugo Rodrigues/Documents/TF_Targets/Promoters/CorkOak_promoters_1kb.bed", sep="\t", names=["Scafold","Start","End","Gene","Strand"])

#Get scafold rows
scafolds = promoter_sequences.loc[::2]
scafolds.reset_index(drop = True, inplace = True)
scafolds["New_Scafolds"] = scafolds["Scafolds"].str[1:]
scafolds = scafolds.drop("Scafolds",axis = 1)

#Get sequences rows
sequences = promoter_sequences.loc[1::2]
sequences.reset_index(drop = True, inplace = True)

#Join them in 2 separate columns
new_promoter_sequences = pd.merge(scafolds, sequences, left_index = True, right_index = True)
new_promoter_sequences.columns = ["Scafold","Sequence"]

#Fix the promoter column and obtain the treated nice_promoter_sequences pandas dataframe
nice_promoter_sequences = pd.concat([new_promoter_sequences, new_promoter_sequences["Scafold"].str.split(':', expand = True)], axis = 1)
nice_promoter_sequences = nice_promoter_sequences.drop(nice_promoter_sequences.columns[[0,3]], axis = 1)
nice_promoter_sequences.columns = ["Sequence", "Scafold"]

#Get scafold-gene link
scafold_gene_link = promoter_co[["Scafold", "Gene"]]
genes = scafold_gene_link["Gene"]
query_promoter_sequences = pd.merge(nice_promoter_sequences, genes, left_index = True, right_index = True)
#query_promoter_sequences.to_csv("C:/Users/Hugo Rodrigues/Documents/TF_Targets/Promoters/Promoter_sequences_df.csv", index=False)
