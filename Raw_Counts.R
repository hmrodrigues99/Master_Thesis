#This script combines the BioProject_Raw_Counts.txt files into a single Raw_Counts dataset.

table1 <- read.delim(file = "~/Coexpression Network/PRJNA690098raw_counts.txt", comment.char = "#", row.names = 1)
table2 <- read.delim(file = "~/Coexpression Network/PRJNA650215raw_counts.txt", comment.char = "#", row.names = 1)
table3 <- read.delim(file = "~/Coexpression Network/PRJ_stems2021raw_counts.txt", comment.char = "#", row.names = 1)
table4 <- read.delim(file = "~/Coexpression Network/PRJNA347903raw_counts.txt", comment.char = "#", row.names = 1)
table5 <- read.delim(file = "~/Coexpression Network/PRJNA392919raw_counts.txt", comment.char = "#", row.names = 1)
table6 <- read.delim(file = "~/Coexpression Network/PRJEB33874raw_counts.txt", comment.char = "#", row.names = 1)

#Cleaning Columns and Columns names
table1_cut <- table1[-c(1:5)]
names(table1_cut) <- gsub(x = names(table1_cut), pattern = "run3.STAR.", replacement = "")
names(table1_cut) <- gsub(x = names(table1_cut), pattern = "run2.STAR.", replacement = "")
names(table1_cut) <- gsub(x = names(table1_cut), pattern = "_R1_merged.filteredAligned.sortedByCoord.out.bam", replacement = "")

table2_cut <- table2[-c(1:5)]
names(table2_cut) <- gsub(x = names(table2_cut), pattern = "aligned_reads.", replacement = "")
names(table2_cut) <- gsub(x = names(table2_cut), pattern = "_pass.filtered.Aligned.out.bam", replacement = "")

table3_cut <- table3[-c(1:5)]                                    
names(table3_cut) <- gsub(x = names(table3_cut), pattern = "X.data.diogolucas.STAR.filtered", replacement = "")
names(table3_cut) <- gsub(x = names(table3_cut), pattern = "Aligned.sortedByCoord.out.bam", replacement = "")

table4_cut <- table4[-c(1:5)]
names(table4_cut) <- gsub(x = names(table4_cut), pattern = "aligned_reads.", replacement = "")
names(table4_cut) <- gsub(x = names(table4_cut), pattern = "_pass.filtered.Aligned.out.bam", replacement = "")

table5_cut <- table5[-c(1:5)]
names(table5_cut) <- gsub(x = names(table5_cut), pattern = "aligned_reads.", replacement = "")
names(table5_cut) <- gsub(x = names(table5_cut), pattern = "_pass.filtered.Aligned.out.bam", replacement = "")

table6_cut <- table6[-c(1:5)]
names(table6_cut) <- gsub(x = names(table6_cut), pattern = "aligned_reads.", replacement = "")
names(table6_cut) <- gsub(x = names(table6_cut), pattern = "_pass.filtered.fq.gz.Aligned.out.bam", replacement = "")

#Merging all tables into the final Raw_Counts dataset
table12 <- transform(merge(table1_cut,table2_cut,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
table123 <- transform(merge(table12,table3_cut,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
table1234 <- transform(merge(table123,table4_cut,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
table12345 <- transform(merge(table1234,table5_cut,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
Raw_Counts <- transform(merge(table12345,table6_cut,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

#Removing the Polen sample (SRR5986741) after interpretation of the DistMatrixHeatmap_Raw_Counts.png because
#it would be a source of unecessary noise.
Raw_Counts$SRR5986741 <- NULL
