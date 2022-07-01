#This script obtains heatmaps with the gene expression of relevant cork network gene modules in the drought project ans seasonal growth project samples

library("DESeq2")
library(tidyverse)
library(ComplexHeatmap)

#Input datasets
Global_Dataset <- read.table("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/Global_Raw_Counts.txt", sep="\t")
Modules_Dataset <- read.table("C:/Users/Hugo Rodrigues/Documents/Module_Analysis/Cork/InfoMap_Clustering_Output/df_LOC_Comunity.txt", header=T)
Cork_Network_Node <- read.table("C:/Users//Hugo Rodrigues/Documents/Cork_network_bundle/DEG_Cork_SLOWP_network04_subset_node.csv", sep = ",", header=T)

#PROCESS DEGS ALL AND DEGS CORK TO JUST HAVE DEGS THAT ARE PRESENT IN MY CORK NETWORK!
DEGs_All <- read.table("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_2/DEG_All.txt", header=T)
Genes_network <- Cork_Network_Node %>% select(name)
DEGs_network <- DEGs_All %>% inner_join(Genes_network, by="name")
DEGs_network <- DEGs_network[!duplicated(DEGs_network$name),]
DEGs_network <- as.data.frame(DEGs_network)
colnames(DEGs_network) <- "name"

DEGs_Cork <- read.table("C:/Users/Hugo Rodrigues/Documents/Cytoscape Networks/DEG_Cork.txt", col.names = "name")
DEGs_Cork_network <- DEGs_Cork %>% inner_join(Genes_network, by="name")
DEGs_Cork_network <- DEGs_Cork_network[!duplicated(DEGs_Cork_network$name),]
DEGs_Cork_network <- as.data.frame(DEGs_Cork_network)
colnames(DEGs_Cork_network) <- "name"

# Generating dummy data to make DeSeq2 datasets.
dummy_meta <- data.frame(N = seq_along(Global_Dataset))
dds <- DESeqDataSetFromMatrix(Global_Dataset, dummy_meta, ~1)

#Pre-filtering - removing genes with low read counts ( < 10 ).
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Variance stabilization to obtain the stabilized data as a matrix.
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- t(assay(vsd))

# Dropping genes that do not vary at all.
vars <- apply(vst, 2, var)
filt_id <- which(is.finite(vars))
vst <- vst[, filt_id]

#Normalized_Dataset
Normalized_Global_Dataset <- t(vst)

# Centering samples around their median to improve reconstruction accuracy. 
#Obtaining Z-scores from normalized values.
means <- apply(Normalized_Global_Dataset, 1, mean)
standarddeviation <- apply(Normalized_Global_Dataset, 1, sd)
Global_Means_Zscore <- sweep(Normalized_Global_Dataset, MARGIN=1, FUN='-', STATS=means)
Global_Means_Zscore <- sweep(Global_Means_Zscore, MARGIN=1, FUN='/', STATS=standarddeviation)
Global_Means_Zscore <- as.data.frame(Global_Means_Zscore)
#write.table(Global_Means_Zscore, file = "C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_1/Global_Means_Zscore")
#write.table(Normalized_Global_Dataset, file = "C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_1/Global_Normalized_Counts")
Global_Means_Zscore$name <- row.names(Global_Means_Zscore)

#Obtain a Gene Module/Cluster composition
Module_1 <- subset(Modules_Dataset,Module=="1")
Module_2 <- subset(Modules_Dataset,Module=="2")
Module_3 <- subset(Modules_Dataset,Module=="3")
Module_4 <- subset(Modules_Dataset,Module=="4")
Module_5 <- subset(Modules_Dataset,Module=="5")

#Filter Clusters only for DEGs_All
cluster_1 <- Module_1 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- cluster_1 %>% inner_join(DEGs_network,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- cluster_2 %>% inner_join(DEGs_network,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- cluster_3 %>% inner_join(DEGs_network,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

cluster_4 <- Module_4 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- cluster_4 %>% inner_join(DEGs_network,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- subset(cluster_4, select = -c(name, Module, Flow))

cluster_5 <- Module_5 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- cluster_5 %>% inner_join(DEGs_network,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- subset(cluster_5, select = -c(name, Module, Flow))

#Filter Clusters only for DEGs_Cork
cluster_1 <- Module_1 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- cluster_1 %>% inner_join(DEGs_Cork_network,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- cluster_2 %>% inner_join(DEGs_Cork_network,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- cluster_3 %>% inner_join(DEGs_Cork_network,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

cluster_4 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- cluster_4 %>% inner_join(DEGs_Cork_network,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- subset(cluster_4, select = -c(name, Module, Flow))

cluster_5 <- Module_5 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- cluster_5 %>% inner_join(DEGs_Cork_network,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- subset(cluster_5, select = -c(name, Module, Flow))

#Change cluster colnames to be compreensible (All samples)
old_names <- colnames(cluster_1)
Project = c(rep("PRJNA690098", 28), rep("PRJNA650215", 11), rep("PRJ_stems2021", 18), rep("PRJNA347903", 4), rep("PRJNA392919", 4), rep("PRJEB33874", 9))
Tissue = c(rep(c("ARoot","ARoot","ARoot","CRoot","CRoot","CRoot"), 3), c("ARoot", "CRoot","ARoot", "ARoot", "ARoot", "CRoot", "CRoot", "CRoot", "ARoot", "CRoot"), rep("Cork", 11), rep(c("Cork", "Phloem", "Xylem"), 6), rep("Embryo", 4), c("Cork","Xylem","Leaf","Innerbark"), rep("Cork", 6), rep("Xylem", 3))
new_names_df <- data.frame(row.names = old_names, Project, Tissue, stringsAsFactors = TRUE)
new_names_df$new_names <- paste(new_names_df$Project,"-", new_names_df$Tissue)
colnames(cluster_1) <- new_names_df$new_names
colnames(cluster_2) <- new_names_df$new_names
colnames(cluster_3) <- new_names_df$new_names
colnames(cluster_4) <- new_names_df$new_names
colnames(cluster_5) <- new_names_df$new_names

#Plot individual Heatmaps for ALL SAMPLES
Heatmap(cluster_1, name = "global_expression", column_title = "Samples", row_title = "cluster_1", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7))
Heatmap(cluster_2, name = "global_expression", column_title = "Samples", row_title = "cluster_2", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7))
Heatmap(cluster_3, name = "global_expression", column_title = "Samples", row_title = "cluster_3", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7))
Heatmap(cluster_4, name = "global_expression", column_title = "Samples", row_title = "cluster_4", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7))
Heatmap(cluster_5, name = "global_expression", column_title = "Samples", row_title = "cluster_5", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7))

#Plot a cluster aggregated Heatmap
All_clusters <- rbind(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
All_clusters_Heatmap <- subset(All_clusters, select = -c(name, Module, Flow))
Heatmap(All_clusters_Heatmap, name = "drought_expression", split = All_clusters$Module, column_title = "Samples", show_row_names = FALSE)