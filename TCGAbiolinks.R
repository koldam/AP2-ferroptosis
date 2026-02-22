library('TCGAbiolinks')
library(dplyr)
library(SummarizedExperiment)
library(GDCRNATools)

query <- GDCquery(project = "TCGA-...",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts")

# download query data using api
GDCdownload(query = query, method = "api")

# prepare data
data <- GDCprepare(query)
data

## generate matrix
counts <- assay(data,1) %>% data.frame()
counts[1:3,1:3]
write.table(counts, file="....txt", sep="\t", col.names = NA)