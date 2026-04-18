#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: get_Pfamscan_lncRNAs.R Pfam_hits lncRNA_npcts_list final_lncRNAs\n")
}

Pfam_hits <- args[1]
lncRNA_npcts_list <- args[2]
final_lncRNAs <- args[3]

data1 <- read.table(Pfam_hits, header=F, sep = '\t')
data2 <- read.table(lncRNA_npcts_list, header=F, sep = '\t')
data3 <- data2[!data2$V1 %in% data1$V1,]
write.table(data3, file = final_lncRNAs , row.names = F, col.names = F, quote = F, sep = '\t')
