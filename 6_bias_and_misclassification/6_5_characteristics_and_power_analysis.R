library(powerSurvEpi)

metadata_ER_POS_SURV     <- readRDS("./output_data/metadata_ER_POS_SURV.RDS")
metadata_ER_POS_REC      <- readRDS("./output_data/metadata_ER_POS_REC.RDS")
refined_data_unique      <- readRDS("./output_data/refined_data_unique.RDS")
metadata_gse_2034_er_pos <- readRDS("./output_data/metadata_gse_2034_er_pos.RDS")
metadata_gse96058_er_pos <- readRDS("./output_data/metadata_gse96058_er_pos.RDS")



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


surv <-
  rbind(
    (
      metadata_ER_POS_SURV %>%
        dplyr::select(EVENT_STAT, EVENT_MON) %>%
        group_by(EVENT_STAT) %>%
        summarize(
          n = n(),
          median_fp = round(median(EVENT_MON), 3),
          min_fp = round(min(EVENT_MON), 3),
         max_fp = round(max(EVENT_MON), 3)
        ) %>%
        mutate(cohort = "METABRIC")
    ),
    
    (
      proof_genes_pt_tcga_survival %>%
        dplyr::select(EVENT_STAT, EVENT_MON) %>%
        group_by(EVENT_STAT) %>%
        summarize(
          n = n(),
          median_fp = median(EVENT_STAT),
          min_fp = min(EVENT_STAT),
          max_fp = max(EVENT_STAT)
        ) %>%
        rename("EVENT_STAT" = "EVENT_STAT") %>%
        mutate(cohort = "TCGA")
    ),
    
    metadata_gse96058_er_pos %>%
      dplyr::select(EVENT_STAT, EVENT_MON) %>%
      group_by(EVENT_STAT) %>%
      summarize(
        n = n(),
        median_fp = round(median(EVENT_MON), 3),
        min_fp = round(min(EVENT_MON), 3),
       max_fp = round(max(EVENT_MON), 3)
      )
    %>%
      mutate(cohort = "GSE96058")
    
    
  ) %>%
  dplyr::rename(
    "N. Patients" = "n",
    "Median follow up" = "median_fp",
    "Minimmum follow up" = "min_fp",
    "Maximmum follow up" = "max_fp",
    "Cohort" = "cohort"
  ) %>% 
  mutate(EVENT_STAT = ifelse(EVENT_STAT == 1, "Deceased", "Alive"))



rec <-
  rbind(
    (
      metadata_ER_POS_REC %>%
        dplyr::select(EVENT_STAT, EVENT_MON) %>%
        group_by(EVENT_STAT) %>%
        summarize(
          n = n(),
          median_fp = round(median(EVENT_MON), 3),
          min_fp = round(min(EVENT_MON), 3),
         max_fp = round(max(EVENT_MON), 3)
        ) %>%
        mutate(cohort = "METABRIC")
    ),
    
    (
      proof_genes_pt_tcga_recurrence %>%
        dplyr::select(EVENT_STAT, EVENT_MON) %>%
        drop_na(EVENT_MON) %>%
        group_by(EVENT_STAT) %>%
        summarize(
          n = n(),
          median_fp = round(median(EVENT_MON), 3),
          min_fp = round(min(EVENT_MON), 3),
          max_fp = round(max(EVENT_MON), 3)
        ) %>%
        rename("EVENT_STAT" = "EVENT_STAT") %>%
        mutate(cohort = "TCGA")
    ),
    
    metadata_gse_2034_er_pos %>%
      dplyr::select(EVENT_STAT, EVENT_MON) %>%
      group_by(EVENT_STAT) %>%
      summarize(
        n = n(),
        median_fp = round(median(EVENT_MON), 3),
        min_fp = round(min(EVENT_MON), 3),
       max_fp = round(max(EVENT_MON), 3)
      )
    %>%
      mutate(cohort = "GSE2034")
    
    
  ) %>%
  dplyr::rename(
    "N. Patients" = "n",
    "Median follow up" = "median_fp",
    "Minimmum follow up" = "min_fp",
    "Maximmum follow up" = "max_fp",
    "Cohort" = "cohort"
  ) %>% 
  mutate(EVENT_STAT = ifelse(EVENT_STAT == 1, "Recurred", "Non recurrent"))



binded <- rbind(
  c(
    "Survival",
    "Survival",
    "Survival",
    "Survival",
    "Survival",
    "Survival"
  ),
  surv,
  c(
    "Recurrence",
    "Recurrence",
    "Recurrence",
    "Recurrence",
    "Recurrence",
    "Recurrence"
  ),
  
  rec
) %>% 
  relocate(Cohort)


binded %>% 
flextable() %>%
  autofit()



power_df <- proof_genes_pt_tcga_recurrence

power_df$group_recoded <- ifelse(power_df$risk_group == "High Risk", "E", "C")

num_E <- sum(proof_genes_pt_tcga_recurrence$risk_group == "High Risk")
num_C <- sum(proof_genes_pt_tcga_recurrence$risk_group == "Low Risk")

power_result <- powerCT(
  formula = Surv(time = EVENT_MON, event = EVENT_STAT) ~ group_recoded,
  dat = power_df,
  nE = num_E,
  nC = num_C,
  RR = 2.0,
  alpha = 0.05
)


print(power_result)

powerEpiCont(
  formula = .pred_linear_pred ~ 1,
  dat = proof_genes_pt_tcga_recurrence,
  var.X1 = ".pred_linear_pred",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)

powerEpiCont(
  formula = .pred_linear_pred ~ 1,
  dat = proof_genes_pt_tcga_survival,
  var.X1 = ".pred_linear_pred",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)



powerEpiCont(
  formula = .pred_linear_pred ~ 1,
  dat = proof_genes_pt_gse96058,
  var.X1 = ".pred_linear_pred",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)


powerEpiCont(
  formula = .pred_linear_pred ~ 1,
  dat = proof_genes_pt_gse2034,
  var.X1 = ".pred_linear_pred",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)


 


powerEpiCont(
  formula = risk_score ~ 1,
  dat = test_data_survival,
  var.X1 = "risk_score",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)


powerEpiCont(
  formula = risk_score ~ 1,
  dat = test_data_recurrence,
  var.X1 = "risk_score",
  var.failureFlag = "EVENT_STAT",
  n = 748,
  theta = exp(1), 
  alpha = 0.05
)
