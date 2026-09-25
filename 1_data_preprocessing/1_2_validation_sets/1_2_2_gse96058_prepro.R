#> In this script we obtain and preprocess the metadata_gse96058 and the counts data

library(GEOquery)
library(dplyr)
library(tidyr)
library(tibble)

# 1.- Download data -----------------------------------------------------------

# 1.1 Download the supplementary file (The actual expression matrix)

directory <- readline("Directory for download of GSE96058 data")

getGEOSuppFiles("GSE96058", baseDir = directory)

# 1.1.2 Read the specific expression file 

raw <- read.csv(paste0(directory, "/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"), 
                row.names = 1, check.names = FALSE)

counts_data_gse96058 <- raw

# 1.2 Download metadata_gse96058

gse <- getGEO("GSE96058", GSEMatrix = TRUE)

# 1.2.2 Asign to object

pheno <- pData(gse[[1]])



# 2.- Preprocess metadata_gse96058 -------------------------------------------------

pheno$characteristics_ch1.3 <- gsub("\\D", "", pheno$characteristics_ch1.3)


metadata_gse96058 <-
  pheno %>%
  mutate(
    tissue = source_name_ch1,
    age = as.numeric(`age at diagnosis:ch1`),
    tumor_size = as.numeric(`tumor size:ch1`),
    lymph_group = `lymph node group:ch1`,
    lymph_status =  `lymph node status:ch1`,
    er_status = as.numeric(`er status:ch1`),
    pgr_status = as.numeric(`pgr status:ch1`),
    her2_status = as.numeric(`her2 status:ch1`),
    ki67_status = as.numeric(`ki67 status:ch1`),
    nhg = as.factor(`nhg:ch1`),
    er_pred_mgc = as.numeric(`er prediction mgc:ch1`),
    # This are predictions made by RNA if mgc its Molecular Gene Classifier which is older than SCN which is Single Sample Classifier (SCAN-B) and if SGC its single gene classifier
    er_pred_sgc = as.numeric(`er prediction sgc:ch1`),
    pgr_pred_mfc = as.numeric(`pgr prediction mgc:ch1`),
    pgr_pred_sgc = as.numeric(`pgr prediction sgc:ch1`),
    her2_pred_mfc = as.numeric(`her2 prediction mgc:ch1`),
    her2_pred_sgc = as.numeric(`her2 prediction sgc:ch1`),
    ki67_pred_mfc = as.numeric(`ki67 prediction mgc:ch1`),
    ki67_pred_sgc = as.numeric(`ki67 prediction sgc:ch1`),
    nhg_pred_mgc = as.numeric(`nhg prediction mgc:ch1`),
    pam50 = as.factor(`pam50 subtype:ch1`),
    os_months = as.numeric(`overall survival days:ch1`) / 30.4166667,
    os_status = as.numeric(`overall survival event:ch1`),
    endocrine_tx = as.numeric(`endocrine treated:ch1`),
    chemo_tx = as.numeric(`chemo treated:ch1`),
    HER2 = her2_pred_sgc, 
    LYMPH = lymph_group,
    PAM50 = pam50,
    AGE = as.numeric(age),
    KI67 = ki67_pred_sgc,
    HORMONE = endocrine_tx,
    CHEMO = chemo_tx
  ) %>%
  dplyr::select(
    -c(
      source_name_ch1,
      characteristics_ch1.2,
      characteristics_ch1.3,
      characteristics_ch1.4,
      characteristics_ch1.5,
      characteristics_ch1.6,
      characteristics_ch1.7,
      characteristics_ch1.8,
      characteristics_ch1.9,
      characteristics_ch1.10,
      characteristics_ch1.11,
      characteristics_ch1.12,
      characteristics_ch1.13,
      characteristics_ch1.14,
      characteristics_ch1.15,
      characteristics_ch1.16,
      characteristics_ch1.17,
      characteristics_ch1.18,
      characteristics_ch1.19,
      characteristics_ch1.20,
      characteristics_ch1.21,
      characteristics_ch1.22,
      characteristics_ch1.23,
      characteristics_ch1.24,
      # Each one of the characteristics_ch1. corresponds to its equivalent in the next lines and both correspond in orther to its characteristic in mutate
      `age at diagnosis:ch1`,
      `tumor size:ch1`,
      `lymph node group:ch1`,
      `lymph node status:ch1`,
      `er status:ch1`,
      `pgr status:ch1`,
      `her2 status:ch1`,
      `ki67 status:ch1`,
      `nhg:ch1`,
      `er prediction mgc:ch1`,
      `er prediction sgc:ch1`,
      `pgr prediction mgc:ch1`,
      `pgr prediction sgc:ch1`,
      `her2 prediction mgc:ch1`,
      `her2 prediction sgc:ch1`,
      `ki67 prediction mgc:ch1`,
      `ki67 prediction sgc:ch1`,
      `nhg prediction mgc:ch1`,
      `pam50 subtype:ch1`,
      `overall survival days:ch1`,
      `overall survival event:ch1`,
      `endocrine treated:ch1`,
      `chemo treated:ch1`
      
    )
  )

counts_data_gse96058[1:5,1:5]


metadata_gse96058_er_pos <-
  metadata_gse96058 %>% 
  filter(er_status == 1) %>% 
  rownames_to_column("id") %>% 
  mutate(EVENT_STAT = os_status,
         EVENT_MON = os_months,
         id = NULL) 

# 3.- Preprocess data -----------------------------------------------------

# 3.1 Match it with metadata_gse96058

# 3.1.1 Identify patients in both sets

common_samples <- intersect(colnames(counts_data_gse96058), metadata_gse96058_er_pos$title)

# 3.1.2 Keep the patients in counts data that also have metadata_gse96058

counts_data_gse96058_erpos <- counts_data_gse96058[, common_samples]

counts_data_gse96058_erpos <- t(counts_data_gse96058_erpos)



saveRDS(counts_data_gse96058_erpos, "./output_data/counts_data_gse96058_erpos.RDS")
saveRDS(metadata_gse96058_er_pos, "./output_data/metadata_gse96058_er_pos.RDS")


rm(list = ls())
gc()
