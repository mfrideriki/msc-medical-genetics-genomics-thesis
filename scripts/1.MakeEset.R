
library(tidyverse)
library(openxlsx)
library(Biobase)

rm(list = ls())
gc()

setwd("./AML_NetBID2/")

table <- read.table("data/original_data/data_mrna_seq_cpm.txt", sep = "\t", header = T)
table <- table[!duplicated(table$Hugo_Symbol),]
table <- subset(table, !is.na(Hugo_Symbol))
rownames(table) <- table$Hugo_Symbol
table <- subset(table, select = -c(Hugo_Symbol, Entrez_Gene_Id))

print(paste("Minimum:", min(table), "Maximum:", max(table)))

metadata <- read.xlsx("data/original_data/data_clinical_sample.xlsx")

metadata <- metadata[metadata$SAMPLE_ID %in% colnames(table),]
rownames(metadata) <- metadata$SAMPLE_ID
table <- table[,colnames(table) %in% metadata$SAMPLE_ID]

table <- table[order(rownames(table)),]
exprs_matrix <- as.matrix(table) 

all(rownames(metadata)==colnames(exprs_matrix)) ## should be TRUE

ExpressionSet <- ExpressionSet(assayData = exprs_matrix, phenoData = new("AnnotatedDataFrame", data = metadata))

saveRDS(ExpressionSet, file = "data/net_eset.rds")
