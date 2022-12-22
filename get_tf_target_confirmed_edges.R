#This script overlaps the interactions in the network where TFs (transcription factors) are involved, and the TF target interactions deposited in a curated database (ConnecTF).
#Inputs
#1) Network node table and edge table from Cytoscape
#2) ConnecTF database network
#3) File with the homology link between the species of interest (in my case Qs.suber) genes and Arabidopsis thaliana genes.
# My work obtained this homology link (quercus_link_thaliana.txt) from a blastp using the quercus_blastp_arabisopsis_link.sh script.

library(dplyr)
library(tidyr)

#RUN THE FOLLOWING LINES ONCE!
#Code used once to get the column with the Arabidopsis IDs into a condensed version (from AT1G123.1 to AT1G123 genes)
#Change the directory here for the Network node table
df_node_at <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Cork_network04_TF_node.csv", header = TRUE, sep=",")
df_node_at <- df_node_at[!apply(df_node_at, 1, function(x) any(x=="")),] 
df_new_node_at <- df_node_at %>% separate(Arabidopsis_gene, c("Arabidopsis_concise", "temp"), remove = FALSE)
df_new_node_at <- subset(df_new_node_at, select = c("name", "Arabidopsis_gene", "Arabidopsis_concise"))
#This will output a .txt file to be imported as a node table into Cytoscape, and add the new Arabidopsis_concise column (use as key column the name column).
write.table(df_new_node_at, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/arabidopsis_concise.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#Getting the Input Network ready to query in ConnectTF -> from "LOC edge LOC" to "AT edge AT" format
#My Cork_network04 loaded into df_Cork_SLOWP_subset has a total of 17,006 edges.
#My resulting df_ConnecTF_network is ready for input in ConnecTF, and has 14,407 edges.

#Network node table from Cytoscape
df1 <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/DEG_Cork_SLOWP/DEG_Cork_SLOWP_network04_subset/DEG_Cork_SLOWP_network04_subset_node.csv", header=TRUE)
nodes_df <- df1 %>% select(name, everything())

#Network edge table from Cytoscape
df2 <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/DEG_Cork_SLOWP/DEG_Cork_SLOWP_network04_subset/DEG_Cork_SLOWP_network04_subset_edge.csv", header=TRUE)
df2$name <- gsub(" \\(interacts with\\) ", ".", df2$name)
edges_df_prov <- df2 %>% separate(name, c('Source', 'Target'))
edges_df <- edges_df_prov %>% select(Source, Target, everything())

df_Cork_SLOWP_subset <- edges_df[ , c('Source','shared.interaction','Target','irp_score')]
names(df_Cork_SLOWP_subset)[names(df_Cork_SLOWP_subset) == "Source"] <- "Quercus_gene1"
names(df_Cork_SLOWP_subset)[names(df_Cork_SLOWP_subset) == "Target"] <- "Quercus_gene2"

#Loading and preprocess the df_quercus_thaliana file which as the LOC link AT
df_quercus_thaliana <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Module_Analysis/new_quercus_linked_thaliana.txt", header=TRUE, sep="\t")

df_quercus_thaliana <- df_quercus_thaliana[!apply(df_quercus_thaliana, 1, function(x) any(x=="")),] 

df_quercus_thaliana <- df_quercus_thaliana %>% separate(hit, c("Arabidopsis_gene_concise1", "Temp"), remove = FALSE)
df_quercus_thaliana <- subset(df_quercus_thaliana, select = -c(Temp, hit))

names(df_quercus_thaliana)[names(df_quercus_thaliana) == "query"] <- "Quercus_gene1"

#Passing the LOC edge LOC network to AT edge AT
df_ConnecTF_input_network <- left_join(df_Cork_SLOWP_subset, df_quercus_thaliana, by="Quercus_gene1")

names(df_quercus_thaliana)[names(df_quercus_thaliana) == "Quercus_gene1"] <- "Quercus_gene2"
names(df_quercus_thaliana)[names(df_quercus_thaliana) == "Arabidopsis_gene_concise1"] <- "Arabidopsis_gene_concise2"

df_ConnecTF_input_network <- left_join(df_ConnecTF_input_network, df_quercus_thaliana, by="Quercus_gene2")

df_ConnecTF_input_network <- subset(df_ConnecTF_input_network, select = c("Arabidopsis_gene_concise1", "shared.interaction", "Arabidopsis_gene_concise2", "irp_score"))

df_ConnecTF_input_network <- na.omit(df_ConnecTF_input_network)

df_ConnecTF_input_network <- df_ConnecTF_input_network[order(-df_ConnecTF_input_network$irp_score), ]

names(df_ConnecTF_input_network)[names(df_ConnecTF_input_network) == "Arabidopsis_gene_concise1"] <- "source"
names(df_ConnecTF_input_network)[names(df_ConnecTF_input_network) == "Arabidopsis_gene_concise2"] <- "target"
names(df_ConnecTF_input_network)[names(df_ConnecTF_input_network) == "shared.interaction"] <- "edge"
df_ConnecTF_input_network["edge"] <- gsub(" ", "_", df_ConnecTF_input_network$edge)

write.table(df_ConnecTF_input_network, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/ConnecTF_input_network.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Ignore for now
df_tf_list <- read.table(file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Ath_TF_list.txt", header = TRUE)
names(df_tf_list)[names(df_tf_list) == "Gene_ID"] <- "Arabidopsis_ID"

df_unique_tf <- left_join(df_new_node_at, df_tf_list, by="Arabidopsis_ID")
df_unique_tf <- na.omit(df_unique_tf)


#Implementing the Output table from ConnecTF (df_cork_tf_target), which has all found TFs and respective targets present in my input network (33,315 edges in 
#df_ConnecTF_network, out of which 24,994 are unique connections, [there are more edges because the same edge can be discovered by different methods]) 
#into the df_ConnecTF_network. The output will be a df_merged, which is the df_ConnecTF_network with the extra columns of "EDGE_TYPE", "Log2FC" and "Pvalue".
#These last 2 columns are rare, because they represent TF target connections that have either up or down expression of the targets.


df_ConnecTF_output <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Cork_allTFs_Targets.csv", header = TRUE, sep=",")
#Did this command for the second ConnecTF run
df_ConnecTF_output <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/ConnecTF_raw_output2.xlsx")
df_ConnecTF_output <- subset(df_cork_tf_target, select = c("gene_id", "EDGE_TYPE", "TARGET", "Log2FC", "Pvalue"))
names(df_ConnecTF_output)[names(df_ConnecTF_output) =="gene_id"] <- "source"
names(df_ConnecTF_output)[names(df_ConnecTF_output) =="TARGET"] <- "target"
df_ConnecTF_output["ConnecTF_TF"] <- rep("yes", dim(df_ConnecTF_output)[[1]])
write.table(df_ConnecTF_output, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/ConnecTF_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#Starting the validation of both occurring edges present in ConnecTF and predicted by co-expression in my Cork network

#Number of edges (TF-Target pairs) found in ConnecTF regardless of discovery method -> duplicate pairs appear because they were found by 2 or more predictive methods
df_check_tf_targets <- subset(df_ConnecTF_output, select = c("source", "target"))
#dim(check_tf_targets)[[1]] -> This gives a total number of 33,315 TF-Target pairs, including duplicates
df_check_tf_targets <- distinct(df_check_tf_targets)
#dim(check_tf_targets)[[1]] -> This gives a total number of 24,994 TF-Target pairs, without duplicates

df_check_cork_network <- subset(df_ConnecTF_input_network, select = c("source", "target"))
df_check_cork_network_inverted <- subset(df_ConnecTF_input_network, select = c("source", "target"))
names(df_check_cork_network_inverted)[names(df_check_cork_network_inverted) =="target"] <- "tmp"
names(df_check_cork_network_inverted)[names(df_check_cork_network_inverted) =="source"] <- "target"
names(df_check_cork_network_inverted)[names(df_check_cork_network_inverted) =="tmp"] <- "source"

#Look at which connections happen simultaneously in the output from ConnecTF (df_check_tf_targets) and already predicted connections in my network (df_check_cork_network).
#Because my network is undirected, and the tf->targets are directional, I matched common edges between tf_targets and my cork_network where genes at the source and target
#switched pales, resulting in the df_check_cork_network and the df_check_cork_network_inverted

df_common_edges <- inner_join(df_check_tf_targets, df_check_cork_network)
df_common_edges_inverted <- inner_join(df_check_tf_targets, df_check_cork_network_inverted)

df_all_common_edges <- rbind(df_common_edges, df_common_edges_inverted)
df_all_common_edges <- distinct(df_all_common_edges)

#The df_merged as 149 AT edge AT connections predicted in my network and confirmed by ConnecTF
df_merged <- inner_join(df_all_common_edges, df_ConnectTF_network, by=c("source", "target"))
dim(df_merged)[[1]]
df_merged <- inner_join(df_merged, df_cork_tf_target, by=c("source", "target"))
dim(df_merged)[[1]]

df_merged <- subset(df_merged, select = -c(edge))
df_merged <- df_merged[ , c("source", "EDGE_TYPE", "target", "irp_score", "Log2FC", "Pvalue")]

write.table(df_merged, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Confirmed_edges_at.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Transform the Confirmed_edges / df_merged from AT genes to LOC genes

names(df_Cork_SLOWP_subset)[names(df_Cork_SLOWP_subset) =="Quercus_gene1"] <- "source"
names(df_Cork_SLOWP_subset)[names(df_Cork_SLOWP_subset) =="Quercus_gene2"] <- "target"

#Step 1 - Change all ATs in the "source" column of df_merged to LOCs
for (i in 1:dim(df_merged)[[1]]) {
  AT_source <- df_merged$source[i]
  LOCs_source <- df_quercus_thaliana[df_quercus_thaliana$Arabidopsis_gene_concise1==AT_source,]
  df_merged_lines <- df_merged[i,]
  df_merged_lines[1,1] <- LOCs_source[1,1] 
  if (dim(LOCs_source)[[1]] > 1) {
    df_merged_lines <- df_merged_lines %>% slice(rep(1:n(), each = dim(LOCs_source)[[1]]))
    df_merged_lines[2,1] <- LOCs_source[2,1]
    if (dim(LOCs_source)[[1]] > 2) {
      df_merged_lines[3,1] <- LOCs_source[3,1]
      if (dim(LOCs_source)[[1]] > 3) {
        df_merged_lines[4,1] <- LOCs_source[4,1]
      }
    }
  }
  if (i == 1) {
    df_cork_confirmed_edges_half <- df_merged_lines
  }
  df_cork_confirmed_edges_half <- rbind(df_cork_confirmed_edges_half, df_merged_lines)
}

#Steo 2 - Change all ATs in the "target" column of df_cork_confirmed_edges_half to LOCs obtaining the df_cork_confirmed_edges
for (i in 1:dim(df_cork_confirmed_edges_half)[[1]]) {
  AT_target <- df_cork_confirmed_edges_half$target[i]
  LOCs_target <- df_quercus_thaliana[df_quercus_thaliana$Arabidopsis_gene_concise1==AT_target,]
  df_merged_lines <- df_cork_confirmed_edges_half[i,]
  df_merged_lines[1,3] <- LOCs_target[1,1]
  if (dim(LOCs_target)[[1]] > 1) {
    df_merged_lines <- df_merged_lines %>% slice(rep(1:n(), each = dim(LOCs_target)[[1]]))
    df_merged_lines[2,3] <- LOCs_target[2,1]
    if (dim(LOCs_target)[[1]] > 2) {
      df_merged_lines[3,3] <- LOCs_target[3,1]
      if (dim(LOCs_target)[[1]] > 3) {
        df_merged_lines[4,3] <- LOCs_target[4,1]
      }
    }
  }
  if (i == 1) {
    df_cork_confirmed_edges <- df_merged_lines
  }
  df_cork_confirmed_edges <- rbind(df_cork_confirmed_edges, df_merged_lines)
}

df_cork_confirmed_edges <- unique(df_cork_confirmed_edges)
df_Cork_SLOWP_subset_TF <- left_join(df_Cork_SLOWP_subset, df_cork_confirmed_edges, by=c("source","target", "irp_score"))
df_Cork_SLOWP_subset_TF <- df_Cork_SLOWP_subset_TF %>% drop_na(EDGE_TYPE)

write.table(df_Cork_SLOWP_subset_TF, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Cork_TF_Confirmed_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Not sure what this does as of now
cytoscape_edge_vec <- c()
for (i in 1:dim(df_Cork_SLOWP_subset_TF)[[1]]) {
  expression <- paste(df_Cork_SLOWP_subset_TF[i,1], " (interacts with) ", df_Cork_SLOWP_subset_TF[i,3], sep = "")
  cytoscape_edge_vec <- append(cytoscape_edge_vec, expression)
}
df_Cork_SLOWP_subset_TF <- cbind(df_Cork_SLOWP_subset_TF, cytoscape_edge_vec)
names(df_Cork_SLOWP_subset_TF)[names(df_Cork_SLOWP_subset_TF) =="EDGE_TYPE"] <- "TF_Target_Pair"

vec1 <- rep("yes", dim(df_Cork_SLOWP_subset_TF)[[1]])
vec2 <- rep("yes", dim(df_Cork_SLOWP_subset_TF)[[1]])
df_Cork_SLOWP_subset_TF["ConnecTF_TF"] <- vec1
df_Cork_SLOWP_subset_TF["ConnecTF_Target"] <- vec2

write.table(df_Cork_SLOWP_subset_TF, file = "C:/Users/Hugo Rodrigues/Documents/TF_Targets/Cork_TF_Confirmed_edges_cyto.txt", sep = "\t", row.names = FALSE, quote = FALSE)

