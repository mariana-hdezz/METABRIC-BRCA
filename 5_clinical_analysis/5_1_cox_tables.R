library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(flextable)
library(officer)

# Load

test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")

train_data_survival <- readRDS("./output_data/train_data_survival.RDS")
train_data_recurrence <- readRDS("./output_data/train_data_recurrence.RDS")

proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")

proof_genes_pt_gse96058  <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")
proof_genes_pt_gse96058_cox  <- readRDS("./output_data/proof_genes_pt_gse96058_cox.RDS")

proof_genes_pt_tcga_cox_survival   <- readRDS("./output_data/proof_genes_pt_tcga_cox_survival.RDS")
proof_genes_pt_tcga_cox_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_cox_recurrence.RDS")

proof_genes_pt_tcga_survival   <- readRDS("./output_data/proof_genes_pt_tcga_survival.RDS")
proof_genes_pt_tcga_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_recurrence.RDS")

proof_genes_pt_cox_survival   <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    test_data <- test_data_recurrence
    analysis <- "recurrence"
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_recurrence
    proof_genes_pt_tcga <- proof_genes_pt_tcga_recurrence
    train_data <- train_data_recurrence
    proof_genes_pt_cox <- proof_genes_pt_cox_recurrence
    
    
  } else{
    test_data <- test_data_survival
    analysis <- "survival"
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_survival
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival
    train_data <- train_data_survival
    proof_genes_pt_cox <- proof_genes_pt_cox_survival
  }
  
  
  # 5.4 Tables
  
  # 5.4.1 Cox regression analysis table with HR
  
  summary_cox <- summary(coxph(surv_obj ~ risk_group, data = test_data))
  
  summary_median <- summary(coxph(surv_obj ~ risk_group_median, data = test_data))
  
  # 6.- C score
  
  # 6.2 Compare concordance
  
  concordance <- concordance(surv_obj ~ risk_score, data = test_data)
  
  c_index_summary <- data.frame(
    C_Index  = concordance$concordance,
    SE       = sqrt(concordance$var),
    data_set = "METABRIC"
  ) %>%
    mutate(
      conf_int_low95  = C_Index - (1.96 * SE),
      conf_int_high95 = C_Index + (1.96 * SE),
      z_stat  = (C_Index - 0.5) / SE,
      p_value = 2 * (1 - pnorm(abs(z_stat)))
    )
  
  clinic_cox <- list(
  summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, data = proof_genes_pt_cox)),
  summary(coxph(surv_obj ~ AGE + LYMPH, data = proof_genes_pt_cox))
  )
  
  # TCGA --------------------------------------------------------------------
  
  # 1.- EXTERNAL VALIDATION ON TCGA
  
  # 2.3.2 Relevel so as to have low risk as reference
  
  proof_genes_pt_tcga$risk_group <- relevel(proof_genes_pt_tcga$risk_group, ref = "Low Risk")
  
  # 2.6 Run Cox again
  
  summary_cox_tcga <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = proof_genes_pt_tcga))
  
  summary_median_tcga <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_tcga))
    
  # 2.7 Calculate the Concordance Index
  
  c_index_results_tcga <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, data = proof_genes_pt_tcga)
  
  
  # 2.8 Table with confidence interval and z stat and estimated p val
  
  c_index_summary_tcga <- data.frame(C_Index = c_index_results_tcga$concordance,
                                     SE = sqrt(c_index_results_tcga$var),
                                     data_set = "TCGA") %>%
    mutate(
      conf_int_low95  = C_Index - (1.96 * SE),
      conf_int_high95 = C_Index + (1.96 * SE),
      z_stat  = (C_Index - 0.5) / SE,
      p_value = 2 * (1 - pnorm(abs(z_stat)))
    )
  
  clinic_cox_tcga <- list(
  summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, data = proof_genes_pt_tcga_cox)),
  summary(coxph(surv_obj ~ AGE + LYMPH  , data = proof_genes_pt_tcga_cox))
  )
  

    if (analysis == "recurrence") {
      c_index_summary_rec <- c_index_summary
      c_index_summary_rec_tcga <- c_index_summary_tcga
      clinic_cox_rec <- clinic_cox
      clinic_cox_rec_tcga <- clinic_cox_tcga
      summary_median_metabric_rec <- summary_median
      summary_median_tcga_rec <- summary_median
    } else{
      c_index_summary_surv <- c_index_summary
      c_index_summary_surv_tcga <- c_index_summary_tcga
      clinic_cox_surv <- clinic_cox
      clinic_cox_surv_tcga <- clinic_cox_tcga
      summary_median_metabric_surv <- summary_median
      summary_median_tcga_surv <- summary_median
    }
  
}


# 1.- GSE2034 ----------------------------------------------------------


# 1.7 Divides time group into less and more than 70 months since if thsi division is not made the proportional hazards assumption is not met

gse2034_split <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ ., 
  data = proof_genes_pt_gse2034,
  cut = 65, 
  episode = "time_group",
  id = "patient_id"
)

# 1.8 Run Cox again

gse2034_split$risk_group <- relevel(gse2034_split$risk_group, ref = "Low Risk")

summary_gse2034 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = gse2034_split))

summary_median_gse2034 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = gse2034_split))


# 1.2 Calculate the Concordance Index

c_index_results_2034 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = proof_genes_pt_gse2034)

# 1.3 Table with confidence interval and z stat and estimated p val

c_index_summary_gse2034 <- data.frame(
  C_Index = c_index_results_2034$concordance,
  SE = sqrt(c_index_results_2034$var),
  data_set = "GSE2034"
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )




# GSE96058 ----------------------------------------------------------------


# 1.4 Cox


proof_genes_pt_gse96058$risk_group <- relevel(proof_genes_pt_gse96058$risk_group, ref = "Low Risk")

summary_gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = proof_genes_pt_gse96058))

summary_median_gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_gse96058))

# 2.5 Calculate the actual Concordance Index

c_index_results_gse96058 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                        data = proof_genes_pt_gse96058)


# 2.6 Table with confidence interval and z stat and estimated p val

c_index_summary_gse96058 <- data.frame(
  C_Index  = c_index_results_gse96058$concordance,
  SE       = sqrt(c_index_results_gse96058$var),
  data_set = "GSE96058"
) %>%
  mutate(
    conf_int_low95  = C_Index - (1.96 * SE),
    conf_int_high95 = C_Index + (1.96 * SE),
    z_stat  = (C_Index - 0.5) / SE,
    p_value = 2 * (1 - pnorm(abs(z_stat)))
  )


clinic_cox_gse96058 <- list(
  summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, data = proof_genes_pt_gse96058_cox)),
  summary(coxph(surv_obj ~ AGE + LYMPH  , data = proof_genes_pt_gse96058_cox))
)



c_index_df_rec <- bind_rows(c_index_summary_rec, c_index_summary_rec_tcga, c_index_summary_gse2034)

c_index_df_rec <- 
  c_index_df_rec %>% 
  relocate(data_set, C_Index, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE)

c_index_df_rec <- 
  c_index_df_rec %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 3)),
    p_value = format(p_value, scientific = TRUE, digits = 3)
  )





c_index_df_surv <- bind_rows(c_index_summary_surv, c_index_summary_surv_tcga, c_index_summary_gse96058)

c_index_df_surv <- 
  c_index_df_surv %>% 
  relocate(data_set, C_Index, conf_int_low95, conf_int_high95, z_stat, p_value) %>% 
  dplyr::select(- SE)

c_index_df_surv <- 
  c_index_df_surv %>% 
  mutate(
    across(where(is.numeric) & !p_value, \(x) round(x, digits = 3)),
    p_value = format(p_value, scientific = TRUE, digits = 3)
  )




summary_median_gse2034_rec <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_gse2034 %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))

summary_median_tcga_rec <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_tcga_recurrence %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))

summary_median_metabric_rec <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = test_data_recurrence %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))



cox_table_rec <- data.frame(
  "Data base" = c("METABRIC", "TCGA", "GSE2034"),
  
  "HR (CI 95%, p-value )" = c(
    paste(
      round(summary_cox$coefficients[2], 3),
      " (", round(summary_cox$conf.int[3], 3),
      " - ", round(summary_cox$conf.int[4], 3),
      ") ",
      ", p-value ", summary_cox$coefficients[5]
    ),
    paste(
      round(summary_cox_tcga$coefficients[2], 3),
      " (", round(summary_cox_tcga$conf.int[3], 3),
      " - ", round(summary_cox_tcga$conf.int[4], 3),
      ") ",
      ", p-value ", summary_cox_tcga$coefficients[5]
    ),
    paste(
      round(summary_gse2034$coefficients[2], 3),
      " (", round(summary_gse2034$conf.int[3], 3),
      " - ", round(summary_gse2034$conf.int[4], 3),
      ") ",
      ", p-value ", summary_gse2034$coefficients[5]
    )
  ),
  
  "C-Index" = c(
    round(summary_cox$concordance[1], 3),
    round(summary_cox_tcga$concordance[1], 3),
    round(summary_gse2034$concordance[1], 3)
  ),
  
  "HR (CI 95%, p-value)" = c(
    paste(
      round(summary_median_metabric_rec$coefficients[2], 3),
      " (", round(summary_median_metabric_rec$conf.int[3], 3),
      " - ", round(summary_median_metabric_rec$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_metabric_rec$coefficients[5]
    ),
    paste(
      round(summary_median_tcga_rec$coefficients[2], 3),
      " (", round(summary_median_tcga_rec$conf.int[3], 3),
      " - ", round(summary_median_tcga_rec$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_tcga_rec$coefficients[5]
    ),
    paste(
      round(summary_median_gse2034$coefficients[2], 3),
      " (", round(summary_median_gse2034$conf.int[3], 3),
      " - ", round(summary_median_gse2034$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_gse2034$coefficients[5]
    )
  ),
  
  "C-Index " = c(
    round(summary_median_metabric_rec$concordance[1], 3),
    round(summary_median_tcga_rec$concordance[1], 3),
    round(summary_median_gse2034$concordance[1], 3)
  )
)



summary_median_gse96058_surv <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_gse96058 %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))

summary_median_tcga_surv <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = proof_genes_pt_tcga_survival %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))

summary_median_metabric_surv <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group_median, data = test_data_survival %>% mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))))



cox_table_surv <- data.frame(
  "Data base" = c("METABRIC", "TCGA", "GSE96058"),
  
  "HR (CI 95%, p-value )" = c(
    paste(
      round(summary_cox$coefficients[2], 3),
      " (", round(summary_cox$conf.int[3], 3),
      " - ", round(summary_cox$conf.int[4], 3),
      ") ",
      ", p-value ", summary_cox$coefficients[5]
    ),
    paste(
      round(summary_cox_tcga$coefficients[2], 3),
      " (", round(summary_cox_tcga$conf.int[3], 3),
      " - ", round(summary_cox_tcga$conf.int[4], 3),
      ") ",
      ", p-value ", summary_cox_tcga$coefficients[5]
    ),
    paste(
      round(summary_gse96058$coefficients[2], 3),
      " (", round(summary_gse96058$conf.int[3], 3),
      " - ", round(summary_gse96058$conf.int[4], 3),
      ") ",
      ", p-value ", summary_gse96058$coefficients[5]
    )
  ),
  
  "C-Index" = c(
    round(summary_cox$concordance[1], 3),
    round(summary_cox_tcga$concordance[1], 3),
    round(summary_gse96058$concordance[1], 3)
  ),
  
  "HR (CI 95%, p-value)" = c(
    paste(
      round(summary_median_metabric_surv$coefficients[2], 3),
      " (", round(summary_median_metabric_surv$conf.int[3], 3),
      " - ", round(summary_median_metabric_surv$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_metabric_surv$coefficients[5]
    ),
    paste(
      round(summary_median_tcga_surv$coefficients[2], 3),
      " (", round(summary_median_tcga_surv$conf.int[3], 3),
      " - ", round(summary_median_tcga_surv$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_tcga_surv$coefficients[5]
    ),
    paste(
      round(summary_median_gse96058$coefficients[2], 3),
      " (", round(summary_median_gse96058$conf.int[3], 3),
      " - ", round(summary_median_gse96058$conf.int[4], 3),
      ") ",
      ", p-value ", summary_median_gse96058$coefficients[5]
    )
  ),
  
  "C-Index " = c(
    round(summary_median_metabric_surv$concordance[1], 3),
    round(summary_median_tcga_surv$concordance[1], 3),
    round(summary_median_gse96058$concordance[1], 3)
  )
)


print(clinic_cox_surv)
print(clinic_cox_surv_tcga)
print(clinic_cox_gse96058)

print(clinic_cox_rec)
print(clinic_cox_rec_tcga)

rbind(
  "Survival",
c_index_df_surv,
"Recurrence",
c_index_df_rec
) %>% 
  dplyr::rename(
    "Cohort" = data_set,
    "C-index" = C_Index,
    "Low 95%" = conf_int_low95,
    "High 95%" = conf_int_high95,
    "Z stat" = z_stat,
    "p value" = p_value
                ) %>% 
  flextable() %>% 
  merge_at(i = 1, j = 1:6, part = "body") %>% 
  merge_at(i = 5, j = 1:6, part = "body") %>% 
  hline(i = 4, border = fp_border(width = 2)) %>% 
  autofit()


rbind(
"Survival",
cox_table_surv,
"Recurrence",
cox_table_rec
) %>% 
  dplyr::rename(
    "Cohort" = Data.base,
    "C - index" = C.Index,
    "HR (CI 95 %, p.value)" = HR..CI.95...p.value..,
    "C - index " = C.Index.,
    "HR (CI 95 %, p.value) " = HR..CI.95...p.value.
  ) %>% 
  flextable() %>% 
  add_header_row(
    values = c("", "Clinic + signature", "Clinic + signature", "Clinic only", "Clinic only"),
    top = TRUE
  ) %>% 
  merge_at(i = 1, j = 2:3, part = "header") %>% 
  merge_at(i = 1, j = 4:5, part = "header") %>%
  merge_at(i = 1, j = 1:5, part = "body") %>% 
  merge_at(i = 5, j = 1:5, part = "body") %>% 
  hline(i = 4, border = fp_border(width = 2)) %>% 
  autofit()

