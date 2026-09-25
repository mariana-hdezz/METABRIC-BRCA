library(timeROC)

# Load

test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
proof_genes_pt_cox_survival <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")

proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")

proof_genes_pt_gse96058  <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")

proof_genes_pt_tcga_survival   <- readRDS("./output_data/proof_genes_pt_tcga_survival.RDS")
proof_genes_pt_tcga_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_recurrence.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    test_data <- test_data_recurrence
    analysis <- "recurrence"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_recurrence
    train_data <- train_data_recurrence
    
    
  } else{
    test_data <- test_data_survival
    analysis <- "survival"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival
    train_data <- train_data_survival
    
  }
  
  
  # 7.- Area under de curve (AUC) per time --------------------------------------------------------
  
  
  # 7.1 Creating object on testing data
  
  time_roc <- timeROC(
    T = test_data$EVENT_MON,
    delta = test_data$EVENT_STAT,
    marker = -test_data$risk_score,
    cause = 1,
    times = c(12, 36, 60, 72, 120),
    # 1y, 3y, 5y, 6y, 10y
    iid = TRUE
  )
  
  
  
  # 2.6.2 Table with confidence interval and z stat and estimated p val
  
  auc_ci <- data.frame(
    AUC  = time_roc$AUC,
    SE   = time_roc$inference$vect_sd_1,
    time = time_roc$times,
    data_set = "METABRIC"
  ) %>%
    mutate(
      conf_int_low95  = AUC - (1.96 * SE),
      conf_int_high95 = AUC + (1.96 * SE),
      z_stat  = (AUC - 0.5) / SE,
      p_value = 2 * (1 - pnorm(abs(z_stat)))
    )
  
  
  # TCGA --------------------------------------------------------------------
  
  
  # 3.1 Area under the curve at 5 time points
  
  res_auc_tcga <- timeROC(
    T = proof_genes_pt_tcga$EVENT_MON,
    delta = proof_genes_pt_tcga$EVENT_STAT,
    marker = -proof_genes_pt_tcga$.pred_linear_pred,
    cause = 1,
    # The EVENT code
    times = c(12, 36, 60, 72, 120),
    # 3, 5, and 10 years
    iid = TRUE
  )
  
  # 3.1.1.2 Table with confidence interval and z stat and estimated p val
  
  auc_ci_tcga <- data.frame(
    AUC  = res_auc_tcga$AUC,
    SE   = res_auc_tcga$inference$vect_sd_1,
    time = res_auc_tcga$times,
    data_set = "TCGA"
  ) %>%
    mutate(
      conf_int_low95  = AUC - (1.96 * SE),
      conf_int_high95 = AUC + (1.96 * SE),
      z_stat  = (AUC - 0.5) / SE,
      p_value = 2 * (1 - pnorm(abs(z_stat)))
    )
  
  if (analysis == "recurrence") {
    auc_ci_rec <- auc_ci
    auc_ci_tcga_rec <- auc_ci_tcga
    saveRDS(time_roc, "./output_data/time_roc_rec.RDS")
    saveRDS(res_auc_tcga, "./output_data/res_auc_tcga_rec.RDS")
    
  } else{
    auc_ci_surv <- auc_ci
    auc_ci_tcga_surv <- auc_ci_tcga
    saveRDS(time_roc, "./output_data/time_roc_surv.RDS") 
    saveRDS(res_auc_tcga, "./output_data/res_auc_tcga_surv.RDS")
    
    }
  
}

# GSE2034 -----------------------------------------------------------------


# 2.- Metric results ------------------------------------------------------

# 2. 1 Area under the curve per time

res_auc_gse2034 <- timeROC(
  T = proof_genes_pt_gse2034$EVENT_MON,
  delta = proof_genes_pt_gse2034$EVENT_STAT,
  marker = -proof_genes_pt_gse2034$.pred_linear_pred,
  cause = 1,
  # The event code
  times = c(12, 36, 60, 120),
  # 3, 5, and 10 years
  iid = TRUE
)

# 2.1.2 Table with confidence interval and z stat and estimated p val

auc_ci_gse2034 <- data.frame(
  AUC  = res_auc_gse2034$AUC,
  SE   = res_auc_gse2034$inference$vect_sd_1,
  time = res_auc_gse2034$times,
  data_set = "GSE2034"
) %>%
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )


# gse96058 ----------------------------------------------------------------


# 2.6 Area under the curve per time

res_auc_gse96058 <- timeROC(
  T = proof_genes_pt_gse96058$EVENT_MON,
  delta = proof_genes_pt_gse96058$EVENT_STAT,
  marker = -proof_genes_pt_gse96058$.pred_linear_pred,
  cause = 1,
  # The event code
  times = c(12, 36, 60, 72),
  # 3, 5, and 6 years
  iid = TRUE
)

# 2.6.2 Table with confidence interval and z stat and estimated p val

auc_ci_gse96058 <- data.frame(
  AUC  = res_auc_gse96058$AUC,
  SE   = res_auc_gse96058$inference$vect_sd_1,
  time = res_auc_gse96058$times,
  data_set = "GSE96058"
) %>%
  mutate(
    conf_int_low95  = AUC - (1.96 * SE),
    conf_int_high95 = AUC + (1.96 * SE),
    z_stat  = (AUC - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )


saveRDS(res_auc_gse2034, "./output_data/res_auc_gse2034.RDS")
saveRDS(res_auc_gse96058, "./output_data/res_auc_gse96058.RDS")

auc_df_rec <- bind_rows(auc_ci_rec, auc_ci_tcga_rec, auc_ci_gse2034)

auc_df_rec <-
  auc_df_rec %>%
  relocate(data_set,
           time,
           AUC,
           conf_int_low95,
           conf_int_high95,
           z_stat,
           p_value) %>%
  dplyr::select(-SE) %>%
  mutate(across(where(is.numeric) &
                  !p_value, \(x) round(x, digits = 3)),
         p_value = format(p_value, scientific = TRUE, digits = 3))



auc_df_surv <- bind_rows(auc_ci_surv, auc_ci_tcga_surv, auc_ci_gse96058)

auc_df_surv <-
  auc_df_surv %>%
  relocate(data_set,
           time,
           AUC,
           conf_int_low95,
           conf_int_high95,
           z_stat,
           p_value) %>%
  dplyr::select(-SE) %>%
  mutate(across(where(is.numeric) &
                  !p_value, \(x) round(x, digits = 3)),
         p_value = format(p_value, scientific = TRUE, digits = 3))

rbind(
"Survival",
auc_df_surv,
"Recurrence",
auc_df_rec
) %>% 
  dplyr::rename(
    "Cohort" = data_set,
    "Time point" = time,
    "Area under the curve" = AUC,
    "Low 95%" = conf_int_low95,
    "High 95%" = conf_int_high95,
    "Z stat" = z_stat,
    "p value" = p_value
  ) %>% 
  flextable() %>% 
  merge_at(i = 1, j = 1:7, part = "body") %>% 
  merge_at(i = 16, j = 1:7, part = "body") %>%
  merge_at(j = 1, i = 2:6, part = "body") %>% 
  hline(i = 6) %>% 
  merge_at(j = 1, i = 7:11, part = "body") %>% 
  hline(i = 11) %>% 
  merge_at(j = 1, i = 12:15, part = "body") %>% 
  hline(i = 15, border = fp_border(width = 2)) %>% 
  merge_at(j = 1, i = 17:21, part = "body") %>% 
  hline(i = 21) %>% 
  merge_at(j = 1, i = 22:26, part = "body") %>% 
  hline(i = 26) %>% 
  merge_at(j = 1, i = 27:30, part = "body") %>% 
  autofit()
