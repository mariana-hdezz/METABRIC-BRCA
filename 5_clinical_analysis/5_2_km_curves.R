library(dplyr)
library(survminer)
library(patchwork)
library(stringr)


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

    
  } else{
    test_data <- test_data_survival
    analysis <- "survival"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival

    
  }

# 5.3 Creating kapan meier curve

# 5.3.1 Based on the Recurrence object compare the risk_groups created previously

fit_km <- survfit(surv_obj ~ risk_group, data = test_data)

# 5.3.2 Plot

km_metabric <- 
  ggsurvplot(fit_km,
           data = test_data,
           pval = TRUE, 
           risk.table = TRUE,
           
           title = paste0(str_to_title(analysis)," ER+ METABRIC"),
           ylab = paste0(str_to_title(analysis), " probability"),
           font.title = 20,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, ifelse(analysis == "recurrence", 200, 300)),
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#c380d3", "#ff89d4"),
)




# TCGA --------------------------------------------------------------------

# 2.4 Fit the KM curve

km_fit_tcga <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = proof_genes_pt_tcga)

# 2.5.2 Plot

km_tcga <- 
  ggsurvplot(km_fit_tcga, 
           data = proof_genes_pt_tcga, 
           pval = TRUE, 
           risk.table = TRUE,
           
           title = paste0("Validation in TCGA (", str_to_title(analysis), ")"),
           font.title = 20,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, ifelse(analysis == "recurrence", 200, 300)),         # Zoom in
           ylim = c(ifelse(analysis == "recurrence", 0.7, 0), 1),
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3,                 # Line size
           palette = c("#E7B800", "#2E9FDF"), # Custom color palette
           
)


if (analysis == "recurrence") {
  
  km_metabric_rec <- km_metabric 
  km_tcga_rec <- km_tcga
  
} else{
  
  km_metabric_surv <- km_metabric
  km_tcga_surv <- km_tcga 
    
}

}

# GSE2034 -----------------------------------------------------------------



# 1.5 Fit the KM curve

km_fit_gse2034 <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = proof_genes_pt_gse2034)

# 1.6 Plot

km_gse2034 <- 
  ggsurvplot(km_fit_gse2034, 
           data = proof_genes_pt_gse2034, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation GSE2034 (Recurrence)",
           font.title = 20,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("High risk", "Low risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 180),         # Zoom in
           ylim = c(0.25, 1),
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3, 
           palette = c("#E41A1C", "#377EB8"),
)


# GSE96058 ----------------------------------------------------------------



# 2.2 Fit the KM curve

km_fit_gse96058 <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = proof_genes_pt_gse96058)

# 2.3 Plot


km_gse96058 <- 
  ggsurvplot(km_fit_gse96058, 
           data = proof_genes_pt_gse96058, 
           pval = TRUE, 
           risk.table = TRUE,
           
           title = "Validation in GSE96058 (Survival)",
           font.title = 20,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("Low risk", "High risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           ylim = c(0.7, 1),
           xlim = c(0, 85),         
           break.time.by = 10,      
           ggtheme = theme_minimal(), 
           
           linewidth = 3,                 
           palette = c("#E41A1C", "#377EB8"), 
           
)

km_plots <- 
  ((((km_metabric_surv$plot | km_tcga_surv$plot) / (km_gse96058$plot + plot_spacer())) /
    ((km_metabric_rec$plot | km_tcga_rec$plot) / (km_gse2034$plot+ plot_spacer()))) + plot_layout(heights = c(1, 1, 3))) 


