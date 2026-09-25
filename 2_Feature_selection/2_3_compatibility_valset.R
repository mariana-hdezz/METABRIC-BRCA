library(dplyr)
library(tidyr)
library(tibble)

# Load signatures

# In 2_2_Boruta.R all signatures where saved but in thsi script we only load the signatures selected for the
# full analysis in the paper

surv_signature <- readRDS("./results/boruta/final_boruta_survival_100_selected_genes.RDS")


recur_signature <- readRDS("./results/boruta/final_boruta_recurrence_500_selected_genes.RDS")

surv_signature <- make.names(surv_signature)
recur_signature <- make.names(recur_signature)

# 1.- Asigning signature genes for GSE2034--------------------------------------------

# load

gene_expres_matrix_gse2034 <- readRDS("./output_data/gene_expres_matrix_gse2034.RDS")

# 1.1 Find genes present in both data sets

colnames(gene_expres_matrix_gse2034) <- make.names(colnames(gene_expres_matrix_gse2034))

common_genes_meta_gse2034 <- intersect(recur_signature, colnames(gene_expres_matrix_gse2034))


# Signature genes in GSE9658 ----------------------------------------------

# Load

counts_data_gse96058_erpos <- readRDS("./output_data/counts_data_gse96058_erpos.RDS")

# 3.2 Find genes present in both data sets

colnames(counts_data_gse96058_erpos) <- make.names(colnames(counts_data_gse96058_erpos))

common_genes_meta_gse96058 <- intersect(surv_signature, colnames(counts_data_gse96058_erpos))



# Common genes with TCGA --------------------------------------------------

# Load

counts_data_tcga <- readRDS("./output_data/counts_data_tcga.RDS")

# 4.- Assigning signature genes --------------------------------------------

# 4.1 Find genes present in both data sets

rownames(counts_data_tcga) <- make.names(rownames(counts_data_tcga))

common_genes_meta_tcga_surv <- intersect(surv_signature, rownames(counts_data_tcga))

common_genes_meta_tcga_rec <- intersect(recur_signature, rownames(counts_data_tcga))



# Finall survival signature

proof_genes_surv <- intersect(common_genes_meta_tcga_surv, common_genes_meta_gse96058)

# Final recurrence signature

proof_genes_rec <- intersect(common_genes_meta_tcga_rec, common_genes_meta_gse2034)

length(surv_signature)

length(proof_genes_surv)

1 - (length(proof_genes_surv) / length(surv_signature))

length(recur_signature)

length(proof_genes_rec)

1 - (length(proof_genes_rec) / length(recur_signature))


saveRDS(proof_genes_surv, "./output_data/proof_genes_surv.RDS")
saveRDS(proof_genes_rec , "./output_data/proof_genes_rec.RDS")

rm(list = ls())
gc
