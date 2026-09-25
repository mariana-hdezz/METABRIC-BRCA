library(dplyr)
library(tidyr)
library(tibble)

# Load

train_data_survival <- readRDS("./output_data/train_data_survival.RDS")
train_data_recurrence <- readRDS("./output_data/train_data_recurrence.RDS")
test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
metadata_ER_POS_REC  <- readRDS("./output_data/metadata_ER_POS_REC.RDS")
metadata_ER_POS_SURV <- readRDS("./output_data/metadata_ER_POS_SURV.RDS")

metadata_gse_2034_er_pos <- readRDS("./output_data/metadata_gse_2034_er_pos.RDS")
proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")

metadata_gse96058_er_pos <- readRDS("./output_data/metadata_gse96058_er_pos.RDS")
proof_genes_pt_gse96058  <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")

refined_data_unique            <- readRDS("./output_data/refined_data_unique.RDS")
proof_genes_pt_tcga_survival   <- readRDS("./output_data/proof_genes_pt_tcga_survival.RDS")
proof_genes_pt_tcga_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_recurrence.RDS")

final_fit_survival <- readRDS("./output_data/final_fit_survival.RDS")
true_cut_survival <- readRDS("./output_data/true_cut_survival.RDS")

final_fit_recurrence <- readRDS("./output_data/final_fit_recurrence.RDS")
true_cut_recurrence <- readRDS("./output_data/true_cut_recurrence.RDS")

proof_genes_surv <- readRDS("./output_data/proof_genes_surv.RDS")
proof_genes_rec <- readRDS("./output_data/proof_genes_rec.RDS")

true_cut_survival <- readRDS("./output_data/true_cut_survival.RDS")
true_cut_recurrence <- readRDS("./output_data/true_cut_recurrence.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    test_data <- test_data_recurrence
    ml_metadata <- metadata_ER_POS_REC
    analysis <- "recurrence"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_recurrence
    train_data <- train_data_recurrence
    proof_genes <- proof_genes_rec
    final_fit <- final_fit_recurrence
    true_cut <- true_cut_recurrence
    
    
  }else{
    test_data <- test_data_survival
    ml_metadata <- metadata_ER_POS_SURV
    analysis <- "survival"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival
    train_data <- train_data_survival
    proof_genes <- proof_genes_surv
    final_fit <- final_fit_survival
    true_cut <- true_cut_survival
    
  }
  
  
  # 5.- Testing -------------------------------------------------------------
  
  # 5.1 Predictions on test data
  
  test_pred <- predict(final_fit, new_data = test_data, type = "linear_pred")
  
  # 5.2 Creating groups for Kaplan-Meier curves
  
  # 5.2.1 Creating column on test data with its prediction
  
  test_data$risk_score <- test_pred$.pred_linear_pred
  
  # 5.2.2 Dividing the groups by median so as to establish a high and low risk and create a column
  
  
  test_data <-
    test_data %>%
    mutate(risk_group =  factor(ifelse(
      risk_score < true_cut$cutpoint[1, 1], "High", "Low"
    )),
    risk_group_median = factor(ifelse(
      risk_score < median(train_data$risk_score), "High", "Low"
    )))
  
  
  
  test_data$risk_group <- relevel(test_data$risk_group, ref = "Low")
  
  
  # 8.2 Create object with parameters to evaluate on the cox model
  
  proof_genes_pt_cox <-
    test_data %>%
    as.data.frame() %>%
    rownames_to_column("PATIENT_ID") %>%
    left_join(ml_metadata, by = "PATIENT_ID", suffix = c("", ".y")) %>% # Join with metadata and eliminate duplicates
    dplyr::select(-ends_with(".y")) %>%
    column_to_rownames("PATIENT_ID") %>%
    mutate(
      SCORE = risk_score # This one is the score of the model
      
    ) %>%
    dplyr::select(
      all_of(proof_genes),
      surv_obj,
      AGE,
      LYMPH,
      HER2,
      MENO,
      HORMONE,
      CHEMO,
      SURGERY,
      PAM50,
      INTCLUST,
      EVENT_STAT,
      EVENT_MON,
      NPI,
      HIST,
      RADIO,
      SCORE
    ) %>%
    na.omit()
  

  saveRDS(
    test_data,
    paste0(
      "./output_data/",
      "test_data_",
      tolower(analysis),
      ".RDS"
    )
  )
  saveRDS(
    proof_genes_pt_cox,
    paste0(
      "./output_data/",
      "proof_genes_pt_cox_",
      tolower(analysis),
      ".RDS"
    )
  )
  
  
  
  
  if (i == 1) {
    # gse2034 -----------------------------------------------------------------

    # 1.1 Predict on GSE2034 set

    gse2034_results <- predict(final_fit, new_data = proof_genes_pt_gse2034, type = "linear_pred") %>%
      bind_cols(proof_genes_pt_gse2034 %>%
                  rownames_to_column("file_name"))

    # 1.4 Create risk groups based on the median of the predictions

    gse2034_results <- gse2034_results %>%
      mutate(risk_group = as.factor(
        ifelse(
          .pred_linear_pred < true_cut$cutpoint$cutpoint[1],
          "High Risk",
          "Low Risk"
        )
      ),
      risk_group_median = factor(ifelse(
        .pred_linear_pred < median(train_data$risk_score), "High", "Low"
      )))

    proof_genes_pt_gse2034 <-
      proof_genes_pt_gse2034 %>%
      rownames_to_column("file_name") %>%
      left_join(gse2034_results, by = "file_name", suffix = c("", ".drop")) %>%
      dplyr::select(-ends_with(".drop"))



    # 3.1 Asigning the correspondant metadata to the tested patients
    proof_genes_pt_gse2034_cox <-
      proof_genes_pt_gse2034 %>%
      as.data.frame() %>%
      left_join(metadata_gse_2034_er_pos,
                by = "file_name",
                suffix = c("", ".y")) %>%
      dplyr::select(-ends_with(".y")) %>%
      column_to_rownames("file_name") %>%
      mutate(SCORE = gse2034_results$.pred_linear_pred) %>%
      dplyr::select(all_of(proof_genes),
                    surv_obj,
                    SCORE,
                    EVENT_MON,
                    EVENT_STAT) %>%
      na.omit()

    saveRDS(
      proof_genes_pt_gse2034,
      paste0(
        "./output_data/",
        "proof_genes_pt_gse2034",
        ".RDS"
      )
    )
    saveRDS(
      proof_genes_pt_gse2034_cox,
      paste0(
        "./output_data/",
        "proof_genes_pt_gse2034_cox",
        ".RDS"
      )
    )


  } else{

    # GSE96058 ----------------------------------------------------------------
    # 1.3 Predict

    gse96058_results <- predict(final_fit, new_data = proof_genes_pt_gse96058, type = "linear_pred") %>%
      bind_cols(proof_genes_pt_gse96058 %>%
                  rownames_to_column("title"))

    # 2.1 Create risk groups based on the median of the predictions

    gse96058_results <- gse96058_results %>%
      mutate(risk_group = as.factor(
        ifelse(
          .pred_linear_pred < true_cut$cutpoint$cutpoint[1],
          "High Risk",
          "Low Risk"
        )
      ),
      risk_group_median = factor(ifelse(
        .pred_linear_pred < median(train_data$risk_score), "High", "Low"
      )))

    proof_genes_pt_gse96058 <-
      proof_genes_pt_gse96058 %>%
      rownames_to_column("title") %>%
      left_join(gse96058_results, by = "title", suffix = c("", ".drop")) %>%
      dplyr::select(-ends_with(".drop"))


    # 3.1 Prepare data for multivariate cox

    proof_genes_pt_gse96058_cox <-
      proof_genes_pt_gse96058 %>%
      as.data.frame() %>%
      left_join(metadata_gse96058_er_pos,
                by = "title",
                suffix = c("", ".y")) %>% # Join with metadata
      dplyr::select(-ends_with(".y")) %>%
      column_to_rownames("title") %>%
      mutate(SCORE = gse96058_results$.pred_linear_pred,
             RISK = gse96058_results$risk_group) %>%
      dplyr::select(
        all_of(proof_genes),
        surv_obj,
        LYMPH,
        PAM50,
        AGE,
        HER2,
        KI67,
        SCORE,
        RISK,
        EVENT_STAT,
        EVENT_MON,
        HORMONE,
        CHEMO
      ) %>%
      na.omit()


    saveRDS(
      proof_genes_pt_gse96058,
      paste0(
        "./output_data/",
        "proof_genes_pt_gse96058",
        ".RDS"
      )
    )
    saveRDS(
      proof_genes_pt_gse96058_cox,
      paste0(
        "./output_data/",
        "proof_genes_pt_gse96058_cox",
        ".RDS"
      )
    )


  }
  

    # TCGA --------------------------------------------------------------------
  
  # 1.1 Predict with the parameters from the fit to the TCGA data
  
  tcga_results <- predict(final_fit, new_data = proof_genes_pt_tcga, type = "linear_pred") %>%
    bind_cols((proof_genes_pt_tcga %>%
                 rownames_to_column("sampleID")))
  # 2.- Divide by risk groups
  
  # 2.3 Create risk groups based on the median of the predictions or the cutpoint
  
  tcga_results <-
    tcga_results %>%
    mutate(risk_group = factor(
      ifelse(
        .pred_linear_pred < true_cut$cutpoint$cutpoint[1],
        "High Risk",
        "Low Risk"
      )
    ),
    risk_group_median = factor(ifelse(
      .pred_linear_pred < median(train_data$risk_score), "High", "Low"
    )))
  
  # 4.1 Prepare objects
  
  # 4.1.1 First merge with score to then merge with the metadata
  
  proof_genes_pt_tcga <-
    proof_genes_pt_tcga %>%
    rownames_to_column("sampleID") %>%
    left_join(tcga_results, by = "sampleID", suffix = c("", ".drop")) %>%
    dplyr::select(-ends_with(".drop"))
  
  # 4.1.2 Merge with metadata and select parameters to evaluate
  
  proof_genes_pt_tcga_cox <-
    proof_genes_pt_tcga %>%
    left_join(refined_data_unique, by = "sampleID") %>%
    as.data.frame() %>%
    column_to_rownames("sampleID") %>%
    mutate(SCORE = tcga_results$.pred_linear_pred, ) %>%
    dplyr::select(
      all_of(proof_genes),
      surv_obj,
      LYMPH,
      PAM50,
      AGE,
      SCORE,
      EVENT_STAT,
      EVENT_MON,
      RADIO,
      NEO,
      OTHER_TX,
      TARG_TX,
      SURGERY,
      risk_group,
      HER2,
      HIST,
      MENO,
      INTCLUST
    )
  
  
  saveRDS(
    proof_genes_pt_tcga,
    paste0(
      "./output_data/",
      "proof_genes_pt_tcga_",
      tolower(analysis),
      ".RDS"
    )
  )
  saveRDS(
    proof_genes_pt_tcga_cox,
    paste0(
      "./output_data/",
      "proof_genes_pt_tcga_cox_",
      tolower(analysis),
      ".RDS"
    )
  )
  
  
}


rm(list = ls())
gc()
