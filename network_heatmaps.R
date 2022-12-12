#This script obtains heatmaps with the gene expression of relevant cork network gene modules in the drought project ans seasonal growth project samples

library("DESeq2")
library(tidyverse)
library(ComplexHeatmap)

#Input datasets
Global_Dataset <- read.table("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/Global_Raw_Counts.txt", sep="\t")
Modules_Dataset <- read.table("C:/Users/Hugo Rodrigues/Documents/Module_Analysis/Cork/InfoMap_Clustering_Output/df_LOC_Comunity.txt", header=T)
Cork_Network_Node <- read.table("C:/Users//Hugo Rodrigues/Documents/Cork_network_bundle/DEG_Cork_SLOWP_network04_subset_node.csv", sep = ",", header=T)
Link_LocAT <- read.table("C:/Users/Hugo Rodrigues/Documents/TF_Targets/arabidopsis_concise.txt", header=T)
TF_Family <- read.table("C:/Users/Hugo Rodrigues/Documents/TF_Targets/Ath_TF_list.txt", header=T)


#PROCESS DEGS ALL AND DEGS CORK TO JUST HAVE DEGS THAT ARE PRESENT IN MY CORK NETWORK!
DEGs_All <- read.table("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_2/DEG_All.txt", header=T)
Genes_network <- Cork_Network_Node %>% select(name)
DEGs_network <- DEGs_All %>% inner_join(Genes_network, by="name")
DEGs_network <- DEGs_network[!duplicated(DEGs_network$name),]
DEGs_network <- as.data.frame(DEGs_network)

DEGs_Cork <- read.table("C:/Users/Hugo Rodrigues/Documents/Cytoscape Networks/DEG_Cork.txt", col.names = "name")
DEGs_Cork_network <- DEGs_Cork %>% inner_join(Genes_network, by="name")
DEGs_Cork_network <- DEGs_Cork_network[!duplicated(DEGs_Cork_network$name),]
DEGs_Cork_network <- as.data.frame(DEGs_Cork_network)

#Sub-Datasets
DStems <- Global_Dataset[ ,40:57]
SStems <- Global_Dataset[ ,29:39]

# Generating dummy data to make DeSeq2 datasets.
dummy_meta <- data.frame(N = seq_along(Global_Dataset))
Ddummy_meta <- data.frame(N = seq_along(DStems))
Sdummy_meta <- data.frame(N = seq_along(SStems))
dds <- DESeqDataSetFromMatrix(Global_Dataset, dummy_meta, ~1)
Ddds <- DESeqDataSetFromMatrix(DStems, Ddummy_meta, ~1)
Sdds <- DESeqDataSetFromMatrix(SStems, Sdummy_meta, ~1)

#Pre-filtering - removing genes with low read counts ( < 10 ).
keep <- rowSums(counts(dds)) >= 10
Dkeep <- rowSums(counts(Ddds)) >= 10
Skeep <- rowSums(counts(Sdds)) >= 10
dds <- dds[keep,]
Ddds <- Ddds[Dkeep,]
Sdds <- Sdds[Skeep,]

# Variance stabilization to obtain the stabilized data as a matrix.
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- t(assay(vsd))
Dvsd <- varianceStabilizingTransformation(Ddds, blind=TRUE)
Dvst <- t(assay(Dvsd))
Svsd <- varianceStabilizingTransformation(Sdds, blind=TRUE)
Svst <- t(assay(Svsd))

# Dropping genes that do not vary at all.
vars <- apply(vst, 2, var)
Dvars <- apply(Dvst, 2, var)
Svars <- apply(Svst, 2, var)

filt_id <- which(is.finite(vars))
Dfilt_id <- which(is.finite(Dvars))
Sfilt_id <- which(is.finite(Svars))

vst <- vst[, filt_id]
Dvst <- Dvst[, Dfilt_id]
Svst <- Svst[, Sfilt_id]

#Normalized Datasets (DStems - stems drought project | SStems - seasonal growth project)
Normalizaded_Global_Dataset <- t(vst)
Normalized_DStems <- t(Dvst)
Normalized_SStems <- t(Svst)

#Calculate replicate means and compact columns in both Sub-Datasets
#Drough Project Sub-Dataset
Xylem_Control <- rowMeans(Normalized_DStems[,c(3,6,9)])
Xylem_Drought <- rowMeans(Normalized_DStems[, c(12,15,18)])
Phloem_Control <- rowMeans(Normalized_DStems[, c(2,5,8)])
Phloem_Drought <- rowMeans(Normalized_DStems[, c(11,14,17)])
Cork_Control <- rowMeans(Normalized_DStems[, c(1,4,7)])
Cork_Drought <- rowMeans(Normalized_DStems[, c(10,13,16)])

Normalized_DStems <- cbind(Normalized_DStems, Xylem_Control, Xylem_Drought, Phloem_Control, Phloem_Drought, Cork_Control, Cork_Drought)
DStems_Means <- Normalized_DStems[ ,19:24]

#Seasonal Growth Project Sub-Dataset
Cork_April <- rowMeans(Normalized_SStems[ ,1:4])
Cork_June <- rowMeans(Normalized_SStems[ ,5:7])
Cork_July <- rowMeans(Normalized_SStems[ ,8:11])

Normalized_SStems <- cbind(Normalized_SStems, Cork_April, Cork_June, Cork_July)
SStems_Means <- Normalized_SStems[ ,12:14]

# Centering samples around their median to improve reconstruction accuracy. 
#Obtaining Z-scores from normalized values.
means <- apply(Normalized_Global_Dataset, 1, mean)
standarddeviation <- apply(Normalized_Global_Dataset, 1, sd)
Global_Means_Zscore <- sweep(Normalized_Global_Dataset, MARGIN=1, FUN='-', STATS=means)
Global_Means_Zscore <- sweep(Global_Means_Zscore, MARGIN=1, FUN='/', STATS=standarddeviation)
Global_Means_Zscore <- as.data.frame(Global_Means_Zscore)

Dmeans <- apply(DStems_Means, 1, mean)
Dstandarddeviation <- apply(DStems_Means, 1, sd)
DStems_Means_Zscore <- sweep(DStems_Means, MARGIN=1, FUN='-', STATS=Dmeans)
DStems_Means_Zscore <- sweep(DStems_Means_Zscore, MARGIN=1, FUN='/', STATS=Dstandarddeviation)
DStems_Means_Zscore <- as.data.frame(DStems_Means_Zscore)

Smeans <- apply(SStems_Means, 1, mean)
Sstandarddeviation <- apply(SStems_Means, 1, sd)
SStems_Means_Zscore <- sweep(SStems_Means, MARGIN=1, FUN='-', STATS=Smeans)
SStems_Means_Zscore <- sweep(SStems_Means_Zscore, MARGIN=1, FUN='/', STATS=Sstandarddeviation)
SStems_Means_Zscore <- as.data.frame(SStems_Means_Zscore)

#Smedians <- apply(Svst, 1, median)
#Svst <- sweep(Svst, MARGIN=1, FUN='-', STATS=Smedians)

Global_Means_Zscore$name <- row.names(Global_Means_Zscore)
#write.table(Global_Means_Zscore, file = "C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_1/Global_Means_Zscore")
DStems_Means_Zscore$name <- row.names(DStems_Means_Zscore)
#write.table(DStems_Means_Zscore, file = "C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_1/Drought_Means_Zscore")
SStems_Means_Zscore$name <- row.names(SStems_Means_Zscore)
#write.table(SStems_Means_Zscore, file = "C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_1/Seasonal_Growth_Means_Zscore")

#Obtain a Gene Module/Cluster composition
Module_1 <- subset(Modules_Dataset,Module=="1")
Module_2 <- subset(Modules_Dataset,Module=="2")
Module_3 <- subset(Modules_Dataset,Module=="3")
Module_4 <- subset(Modules_Dataset,Module=="4")
Module_5 <- subset(Modules_Dataset,Module=="5")

#Obtain Clusters of genes for ALL SAMPLES
cluster_1 <- Module_1 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

cluster_4 <- Module_4 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- subset(cluster_4, select = -c(name, Module, Flow))

cluster_5 <- Module_5 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- subset(cluster_5, select = -c(name, Module, Flow))

#Filter Clusters only for DEGs_All
cluster_1 <- Module_1 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- cluster_1 %>% inner_join(DEGs_All,by="name")
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- cluster_2 %>% inner_join(DEGs_All,by="name")
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- cluster_3 %>% inner_join(DEGs_All,by="name")
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

#Filter Clusters only for DEGs_Cork
cluster_1 <- Module_1 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- cluster_1 %>% inner_join(DEGs_Cork,by="name")
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- cluster_2 %>% inner_join(DEGs_Cork,by="name")
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(Global_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- cluster_3 %>% inner_join(DEGs_Cork,by="name")
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

#Plot individual Heatmaps for ALL SAMPLES
Heatmap(cluster_1, name = "global_expression", column_title = "Samples", row_title = "cluster_1", show_row_names = FALSE, column_names_gp = gpar(fontsize = 6))
Heatmap(cluster_2, name = "global_expression", column_title = "Samples", row_title = "cluster_2", show_row_names = FALSE, column_names_gp = gpar(fontsize = 6))
Heatmap(cluster_3, name = "global_expression", column_title = "Samples", row_title = "cluster_3", show_row_names = FALSE, column_names_gp = gpar(fontsize = 6))
Heatmap(cluster_4, name = "global_expression", column_title = "Samples", row_title = "cluster_4", show_row_names = FALSE, column_names_gp = gpar(fontsize = 6))
Heatmap(cluster_5, name = "global_expression", column_title = "Samples", row_title = "cluster_5", show_row_names = FALSE, column_names_gp = gpar(fontsize = 6))

#Plot a cluster aggregated Heatmap
All_clusters <- rbind(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
All_clusters_Heatmap <- subset(All_clusters, select = -c(name, Module, Flow))
Heatmap(All_clusters_Heatmap, name = "drought_expression", split = All_clusters$Module, column_title = "Samples", show_row_names = FALSE)

#Plot New Heatmaps with New Labels
old_names <- colnames(cluster_1)
Project = c(rep("PRJNA690098", 28), rep("PRJNA650215", 11), rep("PRJ_stems2021", 18), rep("PRJNA347903", 4), rep("PRJNA392919", 4), rep("PRJEB33874", 9))
Tissue = c(rep(c("ARoot","ARoot","ARoot","CRoot","CRoot","CRoot"), 3), c("ARoot", "CRoot","ARoot", "ARoot", "ARoot", "CRoot", "CRoot", "CRoot", "ARoot", "CRoot"), rep("Cork", 11), rep(c("Cork", "Innerbark", "Xylem"), 6), rep("Embryo", 4), c("Cork","Xylem","Leaf","Innerbark"), rep("Cork", 6), rep("Xylem", 3))
new_names_df <- data.frame(row.names = old_names, Project, Tissue, stringsAsFactors = TRUE)
new_names_df$new_names <- paste(new_names_df$Project,"-", new_names_df$Tissue)
colnames(cluster_1) <- new_names_df$new_names
colnames(cluster_2) <- new_names_df$new_names
colnames(cluster_3) <- new_names_df$new_names
colnames(cluster_4) <- new_names_df$new_names
colnames(cluster_5) <- new_names_df$new_names

#Plot and Save new Heatmaps highlighting the first Division
cluster_1 <- as.matrix(cluster_1)
pdf("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/H1_pdf.pdf", width = 9, height = 8, bg = "white", paper = "A4")
Heatmap(cluster_1, name = "Z-Score", column_title = "RNA-seq libraries", row_title = "Module 1", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("Phellem", "Others"), labels_gp = gpar(col = "white", fontize = 10))), column_km = 2)
dev.off()
#ggsave(filename = "H1_gg2", path = "C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/", device = "pdf", width = 185, height = 272, units = "mm", dpi = 300)
cluster_2 <- as.matrix(cluster_2)
pdf("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/H2_pdf.pdf", width = 9, height = 8, bg = "white", paper = "A4")
Heatmap(cluster_2, name = "Z-Score", column_title = "RNA-seq libraries", row_title = "Module 2", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("Others", "Young Root/Embryo"), labels_gp = gpar(col = "white", fontize = 10))), column_km = 2)
dev.off()
cluster_3 <- as.matrix(cluster_3)
pdf("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/H3_pdf.pdf", width = 9, height = 8, bg = "white", paper = "A4")
Heatmap(cluster_3, name = "Z-Score", column_title = "RNA-seq libraries", row_title = "Module 3", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("Young Stem Tissues", "Others"), labels_gp = gpar(col = "white", fontize = 10))), column_km = 2)
dev.off()
cluster_4 <- as.matrix(cluster_4)
pdf("C:/Users/Hugo Rodrigues/Documents/R/Heatmaps/H4_pdf.pdf", width = 9, height = 8, bg = "white", paper = "A4")
Heatmap(cluster_4, name = "Z-Score", column_title = "RNA-seq libraries", row_title = "Module 4", show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4), labels = c("Others", "Roots"), labels_gp = gpar(col = "white", fontize = 10))), column_km = 2)
dev.off()

#Obtain Clusters of genes for the Drought project samples
cluster_1 <- Module_1 %>% inner_join(DStems_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(DStems_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(DStems_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

cluster_4 <- Module_4 %>% inner_join(DStems_Means_Zscore,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- subset(cluster_4, select = -c(name, Module, Flow))

cluster_5 <- Module_5 %>% inner_join(DStems_Means_Zscore,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- subset(cluster_5, select = -c(name, Module, Flow))

#Plot all individual Heatmaps
Heatmap(cluster_1, name = "drought_expression", column_title = "Samples", row_title = "cluster_1", show_row_names = FALSE)
Heatmap(cluster_2, name = "drought_expression", column_title = "Samples", row_title = "cluster_2", show_row_names = FALSE)
Heatmap(cluster_3, name = "drought_expression", column_title = "Samples", row_title = "cluster_3", show_row_names = FALSE)
Heatmap(cluster_4, name = "drought_expression", column_title = "Samples", row_title = "cluster_4", show_row_names = FALSE)
Heatmap(cluster_5, name = "drought_expression", column_title = "Samples", row_title = "cluster_5", show_row_names = FALSE)

#Plot a cluster aggregated Heatmap
All_clusters <- rbind(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
All_clusters_Heatmap <- subset(All_clusters, select = -c(name, Module, Flow))
Heatmap(All_clusters_Heatmap, name = "drought_expression", split = All_clusters$Module, column_title = "Samples", show_row_names = FALSE)

#Obtain Clusters of genes for the Seasonal growth project samples
cluster_1 <- Module_1 %>% inner_join(SStems_Means_Zscore,by="name")
row.names(cluster_1) <- cluster_1$name
cluster_1 <- subset(cluster_1, select = -c(name, Module, Flow))

cluster_2 <- Module_2 %>% inner_join(SStems_Means_Zscore,by="name")
row.names(cluster_2) <- cluster_2$name
cluster_2 <- subset(cluster_2, select = -c(name, Module, Flow))

cluster_3 <- Module_3 %>% inner_join(SStems_Means_Zscore,by="name")
row.names(cluster_3) <- cluster_3$name
cluster_3 <- subset(cluster_3, select = -c(name, Module, Flow))

cluster_4 <- Module_4 %>% inner_join(SStems_Means_Zscore,by="name")
row.names(cluster_4) <- cluster_4$name
cluster_4 <- subset(cluster_4, select = -c(name, Module, Flow))

cluster_5 <- Module_5 %>% inner_join(SStems_Means_Zscore,by="name")
row.names(cluster_5) <- cluster_5$name
cluster_5 <- subset(cluster_5, select = -c(name, Module, Flow))

#Plot all individual Heatmaps
Heatmap(cluster_1, name = "seasonal_expression", column_title = "Samples", row_title = "cluster_1", show_row_names = FALSE)
Heatmap(cluster_2, name = "seasonal_expression", column_title = "Samples", row_title = "cluster_2", show_row_names = FALSE)
Heatmap(cluster_3, name = "seasonal_expression", column_title = "Samples", row_title = "cluster_3", show_row_names = FALSE)
Heatmap(cluster_4, name = "seasonal_expression", column_title = "Samples", row_title = "cluster_4", show_row_names = FALSE)
Heatmap(cluster_5, name = "seasonal_expression", column_title = "Samples", row_title = "cluster_5", show_row_names = FALSE)

#Plot a cluster aggregated Heatmap
All_clusters <- rbind(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
All_clusters_Heatmap <- subset(All_clusters, select = -c(name, Module, Flow))

#Adding color to the aggregated Heatmap
ha = HeatmapAnnotation(df = data.frame(Sample = c("Cork_April", "Cork_June", "Cork_July")))
ha = HeatmapAnnotation(df = data.frame(Sample = c("Cork_April", "Cork_June", "Cork_July")), col = list(Sample = c("Cork_April" = "green", "Cork_June" = "yellow", "Cork_July" = "orange")))

All_clusters_Heatmap <- as.matrix(All_clusters_Heatmap)

#Simpler
Heatmap(All_clusters_Heatmap, name = "seasonal_expression", split = All_clusters$Module, column_title = "Samples", show_row_names = FALSE)

#Adding TF family information to the aggregated Heatmap
colnames(Link_LocAT)[2] <- "TF_ID"
TF_Family_concise <- TF_Family %>% inner_join(Link_LocAT, by="TF_ID")
All_clusters <- All_clusters %>% left_join(TF_Family_concise, by="name")

#Piece of Abstract Art
Heatmap(All_clusters_Heatmap, name = "seasonal_expression", split = All_clusters$Module, top_annotation = ha, show_row_names = FALSE, show_column_names = FALSE) + Heatmap(All_clusters$Module, name = "Module", width = unit(5, "mm"), col = circlize::rand_color(length(unique(All_clusters$Module)))) + Heatmap(All_clusters$Family, name = "TF Family", width = unit(5, "mm"), col = circlize::rand_color(length(unique(All_clusters$Family)))) 
