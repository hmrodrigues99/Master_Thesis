#This script is a simple version of the get_tf_target_confirmed_edges.R

#Necessary libraries
library(dplyr)
library(tidyr)

#Input Cytoscape Node Table here
df_node <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/LOCORK_STUFF/Clean_DEGCork_GCNs/networks_simple/Clean_DEG_GCN4_main_node.csv", header = TRUE, sep=",")
#Input Cytoscape Edge Table here
df_edge <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/LOCORK_STUFF/Clean_DEGCork_GCNs/networks_simple/Clean_DEG_GCN4_main_edge.csv", header = TRUE, sep=",")

df_edge2 <- df_edge[ , c('LOC1', 'interaction', 'LOC2', 'irp_score')]
loc_at <- df_node[ , c('name', 'Arabidopsis_concise')]

names(loc_at)[names(loc_at) == "name"] <- "LOC1"
names(loc_at)[names(loc_at) == "Arabidopsis_concise"] <- "AT1"
df_edge2 <- left_join(df_edge2, loc_at, by="LOC1")

names(loc_at)[names(loc_at) == "LOC1"] <- "LOC2"
names(loc_at)[names(loc_at) == "AT1"] <- "AT2"
df_edge2 <- left_join(df_edge2, loc_at, by="LOC2")

#Remove lines where one or two LOCs dont have a AT match
df_edge2$AT1[df_edge2$AT1 == ""] <- NA
df_edge2$AT2[df_edge2$AT2 == ""] <- NA
df_edge2 <- na.omit(df_edge2)

#Remove duplicate interactions
connectf_query <- df_edge2[ , c('AT1', 'interaction', 'AT2', 'irp_score')]
connectf_query <- connectf_query[!duplicated(connectf_query[c('AT1','interaction','AT2')]),]
write.table(connectf_query, file = "C:/Users/Hugo Rodrigues/Documents/LOCORK_STUFF/Clean_DEGCork_GCNs/networks_simple/Clean_DEG_GCN4_ConnecTF_query.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#   ---   CONNECTF RUN OUTPUT (all_tfs,target_network option) ---

#Resume here:
connectf_output <- read.csv(file = "C:/Users/Hugo Rodrigues/Documents/LOCORK_STUFF/Clean_DEGCork_GCNs/networks_simple/Clean_DEG_GCN4_ConnecTF_output.csv", header = TRUE, sep=",")
connectf_output <- subset(connectf_output, select = c("gene_id", "EDGE_TYPE", "TARGET", "Log2FC", "Pvalue"))

#Merge Source and Target cols into Interaction col
connectf_output <- unite(connectf_output, TF_Target, gene_id,TARGET, sep='_')

#Merge EDGE_TYPE (methods) col into Methods col
connectf_output <- connectf_output %>%
  group_by(TF_Target) %>%
  summarize(EDGE_TYPE=paste0(EDGE_TYPE, collapse = ";"))

#Separate Interaction col into 2 different Source and Target cols
connectf_output <- connectf_output %>%
  separate(TF_Target, c("Source_AT", "Target_AT"), "_")

names(loc_at)[names(loc_at) == "LOC2"] <- "Source_LOC"
names(loc_at)[names(loc_at) == "AT2"] <- "Source_AT"
connectf_output_loc <- left_join(connectf_output, loc_at, by="Source_AT", multiple = 'all')

names(loc_at)[names(loc_at) == "Source_LOC"] <- "Target_LOC"
names(loc_at)[names(loc_at) == "Source_AT"] <- "Target_AT"
connectf_output_loc <- left_join(connectf_output_loc, loc_at, by="Target_AT", multiple = 'all')

#Make LOC interaction
connectf_output_loc$LOC_interaction <- paste(connectf_output_loc$Source_LOC, connectf_output_loc$Target_LOC, sep = ' (interacts with) ')
#And the LOC reverse interaction
connectf_output_loc$LOC_interaction_Inversed <- paste(connectf_output_loc$Target_LOC, connectf_output_loc$Source_LOC, sep = ' (interacts with) ')

#Prepare output table
connectf_output_final <- connectf_output_loc[ , c('LOC_interaction', 'LOC_interaction_Inversed', 'EDGE_TYPE')]
names(connectf_output_final)[names(connectf_output_final) == "EDGE_TYPE"] <- "ConnecTF_Target"

#Write table to Import into Cytoscape!
write.table(connectf_output_final, file = "C:/Users/Hugo Rodrigues/Documents/LOCORK_STUFF/Clean_DEGCork_GCNs/networks_simple/Clean_DEG_GCN4_ConnecTF_Final.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#In the Edge Table Tab, select Import Edge Table, Use Key 'LOC_Interaction' and attribute to add 'ConnecTF_Target' attribute
#Repeat the same Step, with the same table, but in key use 'LOC_Interaction_Inversed' and add 'ConnecTF_Target' atribute


