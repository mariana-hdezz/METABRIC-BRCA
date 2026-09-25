library(survival)
library(dplyr)
library(tidyr)
library(tibble)

# In this file we prepare the data for the cox regression model using onl ER+ patients

# Load

metadata_ER_POS_REC  <- readRDS("./output_data/metadata_ER_POS_REC.RDS")
metadata_ER_POS_SURV <- readRDS("./output_data/metadata_ER_POS_SURV.RDS")
counts_data          <- readRDS("./output_data/counts_data.RDS")

train_rec_id  <- readRDS("./output_data/train_rec_id.RDS")
test_rec_id   <- readRDS("./output_data/test_rec_id.RDS")
train_surv_id <- readRDS("./output_data/train_surv_id.RDS")
test_surv_id  <- readRDS("./output_data/test_surv_id.RDS")

proof_genes_surv <- readRDS("./output_data/proof_genes_surv.RDS")
proof_genes_rec <- readRDS("./output_data/proof_genes_rec.RDS")



for (i in c(1:2)) {
  
  if (i == 1) {
    
    er_patients_metadata <- metadata_ER_POS_REC
    analysis <- "recurrence"
    signature <- proof_genes_rec
    train_id <- train_rec_id
    test_id <- test_rec_id
    
  } else{
    
    er_patients_metadata <- metadata_ER_POS_SURV
    analysis <- "survival"
    signature <- proof_genes_surv
    train_id <- train_surv_id
    test_id <- test_surv_id
    
  }
  
  # 1.- Preparing metadata --------------------------------------------------
  
  ml_metadata <- er_patients_metadata
  
  # 1.2 List of genes to use (check dictionary below to understand the different variables that are used)
  
  proof_genes <- make.names(signature)
  
  
  # 1.3 Object with all the patients ER + and expression of only the genes of interest
  rownames(counts_data) <- make.names(rownames(counts_data))
  proof_genes_pt <- counts_data[proof_genes, er_patients_metadata$PATIENT_ID]
  
  # 1.4  Scaling is done in the cox regression recipe
  
  proof_genes_pt <- t(proof_genes_pt)
  
  # 1.5.1 Check that the patients are in the same order
  
  all(rownames(proof_genes_pt) == er_patients_metadata$PATIENT_ID)
  
  # 1.5.2 Add a column called EVENT_STAT and EVENT_MON to create the surv object
  
  proof_genes_pt <-
    proof_genes_pt %>%
    as.data.frame() %>%
    rownames_to_column("PATIENT_ID") %>%
    left_join(er_patients_metadata, by = "PATIENT_ID") %>%
    column_to_rownames("PATIENT_ID") %>%  # Turn to factor for machine learning
    dplyr::select(all_of(proof_genes), EVENT_MON, EVENT_STAT) %>%
    drop_na() %>%
    mutate(surv_obj = Surv(
      time  = EVENT_MON,
      event = EVENT_STAT,
      type  = "right"
    ))
  
  train_data <-
    proof_genes_pt[rownames(proof_genes_pt) %in% train_id, ]
  
  test_data <-
    proof_genes_pt[rownames(proof_genes_pt) %in% test_id, ]
  
  train_data <- train_data %>%
    mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))
  
  test_data <- test_data %>%
    mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))
  
  proof_genes_pt <- proof_genes_pt %>%
    mutate(across(all_of(proof_genes), ~ as.vector(scale(.x))))
  
  saveRDS(proof_genes_pt, paste0("./output_data/", "proof_genes_pt_" , analysis, ".RDS"))
  saveRDS(train_data,     paste0("./output_data/", "train_data_", analysis, ".RDS"))
  saveRDS(test_data,      paste0("./output_data/", "test_data_", analysis, ".RDS"))
  
}
