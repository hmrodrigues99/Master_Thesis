# Master_Thesis
Collection of custom scripts for the elaboration of my Master Thesis in Bioinformatics and Computational Biology, FCUL 2022, entitled "Identification of gene regulation modules that act in the interaction between cork development and environmental variables"

The code pipeline is as follows:

Task 1) RNA-seq data retrieval and preprocessing:

#RNA-seq libraries retrieval
multi_fastqdump.sh
#Paired reads trimming and quality check
preprocessing.sh
#Create log file of the output from preprocessing.sh
log_reads.sh
#Single reads trimming and quality check
single_preprocessing.sh
#Create log file of the output from single_preprocessing.sh
log_reads_single.sh
#Paired reads alignment to reference genome and obtain read count table using featureCounts
star_script.sh
#Single reads alignment to reference genome and obtain read count table using featureCounts
single_star.sh
#Create log file of the star output of both star_script.sh and single_star_script.sh
log_star.sh
#Create log file of the featureCounts output of both star_script.sh and single_star_script.sh
log_featureCounts.sh
#Merge separate log files into a single one
log_final.sh
#Data integration into a single global Dataset
raw_counts.R
#Data normalization and prepare input into seidr
raw_counts_normalization_seidr.R

Task 2) Prediciton and analysis of the gene co-expression network

#Prediction of the co-expression network
seidr_script.sh
#Isolate the correlation scores of the network and import them into Cytoscape
network_scores.py
#Apply a minimum threshold for the edge scores in the network
network_cut.py
#Plot the distribution of the final correlation scores of the network
network_hist.py
#Obtain a homoly link by blastp between Qs.suber and Ara.thaliana genes
*script_to_do, the commands were simple and done in the shell (prov_name:qs_suber_arabidopsis_thaliana_link.sh)
#Process the blastp results to retain the best hits
blastp_results.py
#Process the network from Cytoscape to input into InfoMap
cytoscape_to_infomap_input.R
#Network clustering using InfoMap
run_infomap.py

#Functional enrichment analysis of infomap network modules using gprofiler
run_gprofiler.R
#Same task as the previous, but applies a initial filter to the modules to only include DEGs (differentialy expressed genes) identified in the Stems Project (Control vs Drought)
filtered_modules_gprofiler.R

#Plot heatmaps with the expression profiles of the top network modules
network_heatmaps.R
#Same task as previous, but applies a initial filter to the modules to only include DEGs identified in the Stems Project
deg_heatmaps.R

Task 3) Identification and analysis of Transcription Factors and putative targets of interest within the network

#Overlaps the ConnecTF database interactions with network predicted interactions to return commonly ocurrying interactions
Get_TF_Target_Confirmed_Edges.R
#Summarize information of network TFs and the genes which they interact with
tf_target_information.py
#Gets the promoter regions of all network genes
*script_to_do, the commands were simple and done in the bash shell (prov_name:get_promoter_regions.sh)
#Perform motif matching between a TF binding site and the promoter region of genes that are connected with that TF in the network (allows input of a TF list)
match_TFBS_promoter.sh

