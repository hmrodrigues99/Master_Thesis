# Master_Thesis
#### Collection of custom scripts for the elaboration of my Master Thesis in Bioinformatics and Computational Biology, FCUL 2022, entitled "Identification of gene regulatory modules acting in the interaction between cork development and environmental variables"

![image](/Thesis_Data_Workflow_Final.png)

The code pipeline is as follows:

### Task 1) RNA-seq data retrieval and preprocessing

**multi_fastqdump.sh** #RNA-seq libraries retrieval

**preprocessing.sh** #Paired reads trimming and quality check

**log_reads.sh** #Create log file of the output from preprocessing.sh

**single_preprocessing.sh** #Single reads trimming and quality check

**log_reads_single.sh** #Create log file of the output from single_preprocessing.sh

**star_script.sh** #Paired reads alignment to reference genome and obtain read count table using featureCounts

**single_star.sh** #Single reads alignment to reference genome and obtain read count table using featureCounts

**log_star.sh** #Create log file of the star output of both star_script.sh and single_star_script.sh

**log_featureCounts.sh** #Create log file of the featureCounts output of both star_script.sh and single_star_script.sh

**log_final.sh** #Merge separate log files into a single one

**Raw_Counts.R** #Data integration into a single global Dataset

**Raw_Counts_normalization.R** (Incorporated into seidr_script.sh) #Data normalization and prepare input into seidr

### Task 2) Prediciton and analysis of the gene co-expression network

**seidr_script.sh** #Prediction of the co-expression network

**network_scores.py** #Isolate the correlation scores of the network and import them into Cytoscape

**network_cut.py** #Apply a minimum threshold for the edge scores in the network

**network_hist.py** #Plot the distribution of the final correlation scores of the network

**quercus_blastp_arabidopsis_link.sh** #Obtain a sequence homology link by blastp between Quercus suber and Arabidopsis thaliana genes

**blastp_results.py** #Process the blastp results to retain the best hits

**cytoscape_to_infomap_input.R** #Process the network from Cytoscape to input into InfoMap

**run_infomap.py** #Network clustering using InfoMap

**run_gprofiler.R** #Functional enrichment analysis of infomap network modules using gprofiler

**filtered_modules_gprofiler.R** #Same task as the previous, but applies a initial filter to the modules to only include DEGs (differentialy expressed genes) identified in the Stems Project (Control vs Drought)

**network_heatmaps.R** #Plot heatmaps with the expression profiles of the top network modules

**deg_heatmaps.R** #Same task as previous, but applies a initial filter to the modules to only include DEGs identified in the Stems Project

### Task 3) Identification and analysis of Transcription Factors and putative targets of interest within the network

**simple_ConnecTF.R** #Overlaps the ConnecTF database interactions with network predicted interactions to return commonly ocurrying interactions

**tf_target_information.py** #Summarize information of network TFs and the genes which they interact with

**get_All_promoter_sequences.sh** #Gets the promoter regions of all cork oak genes

**get_Promoter_sequence_df.py** #Makes a table with scaffold, gene and promoter region columns of cork oak genes present in the network

**match_TFBS_promoter.sh** #Perform motif matching between a TF binding site and the promoter region of genes that are connected with that TF in the network. It incorporates **match_TFBS_promoter.py** to process input files

