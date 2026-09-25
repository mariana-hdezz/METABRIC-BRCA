# In this file we run the external validation for TCGA

library(coxphf)
library(survminer)
library(timeROC)

proof_genes_pt_tcga_survival   <- readRDS("./output_data/proof_genes_pt_tcga_survival.RDS")
proof_genes_pt_tcga_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_recurrence.RDS")
proof_genes_pt_gse2034   <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")
proof_genes_pt_gse96058   <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")

proof_genes_pt_gse2034_cox   <- readRDS("./output_data/proof_genes_pt_gse2034_cox.RDS")
proof_genes_pt_gse96058_cox   <- readRDS("./output_data/proof_genes_pt_gse96058_cox.RDS")
proof_genes_pt_tcga_cox_survival   <- readRDS("./output_data/proof_genes_pt_tcga_cox_survival.RDS")
proof_genes_pt_tcga_cox_recurrence   <- readRDS("./output_data/proof_genes_pt_tcga_cox_recurrence.RDS")

cox_model_gse96058 <- readRDS("./output_data/cox_model_gse96058.RDS")
cox_model_gse2034  <- readRDS("./output_data/cox_model_gse2034.RDS")
cox_model_tcga_rec     <- readRDS("./output_data/cox_model_tcga_rec.RDS")
cox_model_tcga_surv     <- readRDS("./output_data/cox_model_tcga_surv.RDS")

final_fit_survival     <- readRDS("./output_data/final_fit_survival.RDS")
final_fit_recurrence     <- readRDS("./output_data/final_fit_recurrence.RDS")

# 7.- Other scores --------------------------------------------------------

for (i in c(1:2)) {
  if (i == 1) {
    proof_genes_pt_tcga <- proof_genes_pt_tcga_recurrence
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_recurrence
    cox_model_tcga <- cox_model_tcga_rec
    final_fit <- final_fit_recurrence
  } else{
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival
    cox_model_tcga <- cox_model_tcga_surv
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_survival
    final_fit <- final_fit_survival
  }
  
  
  # 7.1 Brier score
  
  eval_results_tcga <- final_fit %>%
    augment(new_data = proof_genes_pt_tcga, eval_time = c(36, 60, 120))
  
  performance_tcga <- eval_results_tcga %>%
    brier_survival(truth = surv_obj, .pred)
  
  print(performance_tcga)
  
  # 7.2 Martingale and Schofeild residuals
  
  cox.zph(cox_model_tcga)
  
  ggcoxzph(cox.zph(cox_model_tcga))
  
  ggcoxdiagnostics(
    cox_model_tcga,
    type = "martingale",
    linear.predictions = FALSE,
    ggtheme = theme_bw()
  )
  
  
  if (i == 1) {
    # 5.- Other scores --------------------------------------------------------
    
    # 5.1 Brier score
    
    eval_results_gse2034 <- final_fit %>%
      augment(new_data = proof_genes_pt_gse2034, eval_time = c(36, 60, 120))
    
    performance_gse2034 <- eval_results_gse2034 %>%
      brier_survival(truth = surv_obj, .pred)
    
    print(performance_gse2034)
    
    # 5.2 Martingale and Schofeild residuals
    
    cox.zph(cox_model_gse2034)
    
    ggcoxzph(cox.zph(cox_model_gse2034))
    
    ggcoxdiagnostics(
      cox_model_gse2034,
      type = "martingale",
      linear.predictions = FALSE,
      ggtheme = theme_bw()
    )
    
    
  } else{
    # 6.- Other evaluations ---------------------------------------------------
    
    # 6.1 Brier score
    
    eval_results_gse96058 <- final_fit %>%
      augment(new_data = proof_genes_pt_gse96058, eval_time = c(36, 60, 120))
    
    # 6.2 Calculate Brier Score
    
    performance_gse96058 <- eval_results_gse96058 %>%
      brier_survival(truth = surv_obj, .pred)
    
    print(performance_gse96058)
    
    # 6.3 Martingale and Schofeild residuals
    
    cox.zph(cox_model_gse96058)
    
    ggcoxzph(cox.zph(cox_model_gse96058))
    
    ggcoxdiagnostics(
      cox_model_gse96058,
      type = "martingale",
      linear.predictions = FALSE,
      ggtheme = theme_bw()
    )
  }
  
}
