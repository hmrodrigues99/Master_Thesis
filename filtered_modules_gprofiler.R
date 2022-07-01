#Run g:Profiler to perform functional enrichment analysis on the top 5 modules filtered for th DEGs occurying in Cork tissue
library(gprofiler2)

outdir <- "C:/Users/Hugo Rodrigues/Documents/Module_Analysis/"

#Reading the output dataframe of InfoMap and the dataframe with the homology link between Quercus LOC and Arabidopsis AT genes
df_LOC_Comunity <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Module_Analysis/Cork/InfoMap_Clustering_Output/df_LOC_Comunity.txt", header=TRUE, sep=" ", comment.char = "#")
df_quercus_thaliana <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Module_Analysis/new_quercus_linked_thaliana.txt", header=TRUE, sep="\t")
#Reading the Cork network (cutoff at 0.4)
Cork_Network_Node <- read.table("C:/Users//Hugo Rodrigues/Documents/Cork_network_bundle/DEG_Cork_SLOWP_network04_subset_node.csv", sep = ",", header=T)
Genes_network <- Cork_Network_Node %>% select(name)
#Reading the list of previously identified differentialy expressed genes in Cork tissue
DEGs_Cork <- read.table("C:/Users/Hugo Rodrigues/Documents/Cytoscape Networks/DEG_Cork.txt", col.names = "name")
DEGs_Cork_network <- DEGs_Cork %>% inner_join(Genes_network, by="name")
DEGs_Cork_network <- DEGs_Cork_network[!duplicated(DEGs_Cork_network$name),]
DEGs_Cork_network <- as.data.frame(DEGs_Cork_network)
colnames(DEGs_Cork_network) <- "Quercus_gene"
#For some weird reason, removing the dot from the genes leads to an error in the g:Profiller run -_("-")_-
#df_quercus_thaliana$Arabidopsis_gene <- gsub("\\.*", "", df_quercus_thaliana$Arabidopsis_gene)
names(df_LOC_Comunity)[names(df_LOC_Comunity) == "name"] <- "Quercus_gene"
names(df_quercus_thaliana)[names(df_quercus_thaliana) == "query"] <- "Quercus_gene"
names(df_quercus_thaliana)[names(df_quercus_thaliana) == "hit"] <- "Arabidopsis_gene"
df_Comunities <- merge(df_LOC_Comunity, df_quercus_thaliana, by="Quercus_gene")
df_Comunities <- df_Comunities[order(df_Comunities$Module),]

nmodules <- tail(df_Comunities$Module, n=1)
loc_genes_per_module <- c()
modules_list <- list()
for(i in 1:5) {
  #get a separate gene list (without duplicates) of each module present in df_CLU (module1, module2, .., modulen) and append it into modules_list
  module <- paste("module", i, sep="_")
  loc_gene_list <- df_LOC_Comunity[which(df_LOC_Comunity$Module == i), ]
  loc_gene_list <- loc_gene_list %>% inner_join(DEGs_Cork_network,by="Quercus_gene")
  loc_genes_per_module <- append(loc_genes_per_module, dim(loc_gene_list)[1])
  gene_list <- df_Comunities[which(df_Comunities$Module == i), ]
  gene_list <- subset(gene_list, select = "Arabidopsis_gene")
  assign(module, unique(gene_list$Arabidopsis_gene))
  tmp <- get(module)
  #Define a minimum gene length of the modules to analyse (default will be 0 / all modules)
  if (length(tmp) > 20) {
    modules_list <- append(modules_list, module)
  }
}

loc_genes_per_module <- loc_genes_per_module[1:length(modules_list)]
print(length(modules_list))
module_names <- unlist(modules_list)
names(modules_list) <- module_names
gem_modules_list <- list()

#Perform functional profilling using g:Profiler with a loop of single module inputs (that meet the minimum gene treshold criteria)
for (i in 1:length(modules_list)) {
  gostres_module <- paste("gostres", names(modules_list)[i], sep="_")
  tmp <- get(modules_list[[i]])
  #gprofiler single query giving all collected go_terms in a single module
  single_module <- gost(query = tmp, organism = "athaliana", ordered_query = FALSE, multi_query = FALSE, significant = FALSE)
  #colecting the enriched go_terms of a single module
  enriched_single_module <- subset(single_module$result, significant=="TRUE")
  #checking if the module is enriched in at least 1 go_term to proceed with the analysis
  if (nrow(enriched_single_module) >= 1) {
    gem_modules_list <- append(gem_modules_list, gostres_module)
    assign(gostres_module, enriched_single_module)
    gostres_gem_module <- enriched_single_module[,c("term_id", "term_name", "p_value", "intersection_size")]
    colnames(gostres_gem_module) <- c("GO.ID", "Description", "p.Val", "Genes")
    gostres_gem_module$FDR <- gostres_gem_module$p.Val
    gostres_gem_module$Phenotype <- "+1"
    gostres_gem_module <- gostres_gem_module[, c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    out_gem <- paste(outdir, "gostres_gem_", names(modules_list)[i], ".txt", sep="")
    write.table(gostres_gem_module, file = out_gem, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

iterations <- length(gem_modules_list)
variables <- 9
output <- matrix(ncol = variables, nrow = iterations)

for(i in 1:iterations ) {
  bp_go_terms_module <- paste("bp_go_terms", gem_modules_list[[i]], sep="_")
  tmp1 <- get(modules_list[[i]])
  tmp2 <- get(gem_modules_list[[i]])
  module_name <- gem_modules_list[[i]]
  module_number <- gsub(".*_.*_", "", module_name)
  at_gene_number <- length(tmp1)
  enriched_GOTerms <- dim(tmp2)[1]
  bp_go_terms <- subset(tmp2, source == "GO:BP")
  bp_go_terms_list <- bp_go_terms$term_id
  out_bp <- paste(outdir, bp_go_terms_module, ".txt", sep="")
  write(bp_go_terms_list, file = out_bp)
  go_terms <- table(tmp2$source)
  go_mf <- 0
  check <- "GO:MF" %in% tmp2$source
  if (check == "TRUE") {
    go_mf <- go_terms[["GO:MF"]]
  }
  go_cc <- 0
  check <- "GO:CC" %in% tmp2$source
  if (check == "TRUE") {
    go_cc <- go_terms[["GO:CC"]]
  }
  go_bp <- 0
  check <- "GO:BP" %in% tmp2$source
  if (check == "TRUE") {
    go_bp <- go_terms[["GO:BP"]]
  }
  kegg <- 0
  check <- "KEGG" %in% tmp2$source
  if (check == "TRUE") {
    kegg <- go_terms[["KEGG"]]
  }
  mirna <- 0
  check <- "MIRNA" %in% tmp2$source
  if (check == "TRUE") {
    mirna <- go_terms[["MIRNA"]]
  }
  wp <- 0
  check <- "WP" %in% tmp2$source
  if (check == "TRUE") {
    wp <- go_terms[["WP"]]
  }
  module_info <- c(module_number, enriched_GOTerms, go_mf, go_cc, go_bp, kegg, mirna, wp, at_gene_number)
  output[i,] <- module_info
}

df_enrichment_data <- data.frame(output)
df_enrichment_data <- setNames(df_enrichment_data,c("Module", "Enriched_GO:Terms", "GO:MF", "GO:CC", "GO:BP", "KEGG", "MIRNA", "WP","At_Gene_Number"))
df_enrichment_data$Loc_Gene_Number <- loc_genes_per_module[1:iterations]
df_enrichment_data <- df_enrichment_data[, c("Module", "Loc_Gene_Number", "At_Gene_Number", "Enriched_GO:Terms", "GO:MF", "GO:CC", "GO:BP", "KEGG", "MIRNA", "WP")]
out_enrichment_data <- paste(outdir, "enrichment_data.csv", sep="")
write.csv(df_enrichment_data, file = out_enrichment_data, row.names = FALSE)
print(df_enrichment_data)
