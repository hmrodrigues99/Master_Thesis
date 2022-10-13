#Run g:Profiler to perform functional enrichment analysis on the top X amount of previously obtained modules.
library(gprofiler2)

outdir <- "C:/Users/Hugo Rodrigues/Documents/Module_Analysis/"

#Reading the output dataframe of InfoMap and the dataframe with the homology link between Quercus LOC and Arabidopsis AT genes
df_LOC_Comunity <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Module_Analysis/Cork/InfoMap_Clustering_Output/df_LOC_Comunity.txt", header=TRUE, sep=" ", comment.char = "#")
df_quercus_thaliana <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Module_Analysis/quercus_arabidopsis_linked.txt", header=TRUE, sep="\t")
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
for(i in 1:nmodules) {
  #get a separate gene list (without duplicates) of each module present in df_CLU (module1, module2, .., modulen) and append it into modules_list
  module <- paste("module", i, sep="_")
  loc_gene_list <- df_LOC_Comunity[which(df_LOC_Comunity$Module == i), ]
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



#Getting Pictures and detailed info on each of the query modules (top_5_modules)
p1 = gostplot(gostres_module_1, interactive = TRUE)
p2 = gostplot(gostres_module_2, interactive = TRUE)
p3 = gostplot(gostres_module_3, interactive = TRUE)
#Saving Images after selecting cool GO:Terms
publish_gostplot(p1, highlight_terms = c("GO:0003824", "GO:0106310", "GO:0031224", "GO:0006629", "GO:0055114", "KEGG:01110"))
publish_gostplot(p2, highlight_terms = c("GO:0071554", "GO:0045229", "GO:0071555", "KEGG:04075"))
publish_gostplot(p3, highlight_terms = c("GO:0009579", "GO:0009534", "GO:0031976", "GO:0042651", "GO:0034357"))
publish_gostplot(p3, highlight_terms = c("GO:0015979", "GO:0019684", "GO:0006091", "KEGG:00195"))

#Perform functional profiling using g:Profiler with multiple queries (top 5 modules of the network)
gostres_multi <- gost(query = top_5_modules , organism = "athaliana", ordered_query = FALSE, multi_query = FALSE)
head(gostres_multi$result)

#Visualizing the comparison between the top 5 modules
p_top5 = gostplot(gostres_multi, interactive = FALSE)
publish_gostplot(p_top5)

#Multiple queries (top 5 modules of the network) and compare the results using the multi_query option
#gostres_multicomp <- gost(query = top_5_modules , organism = "athaliana", ordered_query = FALSE, multi_query = TRUE)
#head(gostres_multicomp$result)
