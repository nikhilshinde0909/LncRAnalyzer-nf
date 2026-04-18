#!/usr/bin/env Rscript

library(dplyr)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: get_CPAT_results.R <CPAT_output> <cutoff_value> <lncRNA_list_out> <NPCTs_list_out>")
}

CPAT_output  <- args[1]
cutoff_value <- args[2]
lncRNA_list  <- args[3]
NPCTs_list   <- args[4]

data <- read.table(CPAT_output, header = TRUE, sep = '\t')

data <- data %>% rownames_to_column("Transcript")

data$coding_prob <- as.numeric(data$coding_prob)

data1 <- read.table(cutoff_value, header = FALSE, sep = '\t')
cutoff <- as.numeric(data1$V1[1])

Lnc <- data[data$coding_prob < cutoff, ]
write.table(Lnc$Transcript, lncRNA_list, row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

NPCTs <- data[data$coding_prob >= cutoff, ]
write.table(NPCTs$Transcript, NPCTs_list, row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
