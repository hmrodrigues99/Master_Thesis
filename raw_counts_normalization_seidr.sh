#!/usr/bin/Rscript

library("argparse")
parser <- ArgumentParser(description="Get the expression.tsv and genes.txt from a Raw_Counts table")
parser$add_argument("-f", "--rawcountsfile", help="Raw_Counts.txt table delimited by spaces")
parser$add_argument("-o", "--output", metavar="STR", help="Output files (default = seidr_output)")
#parser$add_argument("-o1", "--out_expr", metavar="STR", type="character", default='expression.tsv', help="specified name for the expression output file (default = expression.tsv")
#parser$add_argument("-o2", "--out_genes", metavar="STR", type="character", default='genes.txt', help="specified name for the genes list output file (default = genes.txt)")

args <- parser$parse_args()
rawcountsfile <- args$rawcountsfile
output <- args$output

suppressPackageStartupMessages(library("DESeq2"))

Raw_Counts <- read.delim(file = rawcountsfile, comment.char = "#", row.names = 1)

# In order to make a DESeq2 data set, we need some metadata. For now we'll just use some dummy data
dummy_meta <- data.frame(N = seq_along(Raw_Counts))
dds <- DESeqDataSetFromMatrix(Raw_Counts, dummy_meta, ~1)

#Pre-filtering - removing genes with low read counts ( < 10 )
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Now we run variance stabilization and get the stabilized data as a matrix.
# If you have good metadata, you can use the experimental design in the DESeqDataSet and set blind=FALSE here
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- t(assay(vsd))

# Genes that do not vary at all create problems down the line, so it's better to drop them
vars <- apply(vst, 2, var)
filt_id <- which(is.finite(vars))
vst <- vst[, filt_id]

# Let's also center samples around their median, which has been shown to improve reconstruction accuracy
medians <- apply(vst, 1, median)
vst <- sweep(vst, MARGIN=1, FUN='-', STATS=medians)

# MASS's write.matrix function is a bit faster and better suited for our task when compared
# to write.table. Don't forget to unname(), otherwise you will have column headers in the output
#MASS::write.matrix(x=unname(vst), sep='\t', file='expression.tsv')
MASS::write.matrix(x=unname(vst), sep='\t', file=sprintf('./%s/expression.tsv', output))

# Finally, let's write the column headers (== gene names) as a text file
#write(colnames(vst), file='genes.txt')
write(colnames(vst), file=sprintf('./%s/genes.txt', output))
