library(timeROC)
library(dplyr)



test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
proof_genes_pt_cox_survival <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")

proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")

proof_genes_pt_gse96058  <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")

proof_genes_pt_tcga_survival   <- readRDS("./output_data/proof_genes_pt_tcga_survival.RDS")
proof_genes_pt_tcga_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_recurrence.RDS")

res_auc_gse2034  <- readRDS("./output_data/res_auc_gse2034.RDS")
res_auc_gse96058 <- readRDS("./output_data/res_auc_gse96058.RDS")
time_roc_surv     <- readRDS("./output_data/time_roc_surv.RDS")
res_auc_tcga_surv <- readRDS("./output_data/res_auc_tcga_surv.RDS")
time_roc_rec     <- readRDS("./output_data/time_roc_rec.RDS")
res_auc_tcga_rec <- readRDS("./output_data/res_auc_tcga_rec.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    test_data <- test_data_recurrence
    analysis <- "recurrence"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_recurrence
    time_roc <- time_roc_rec
    res_auc_tcga <- res_auc_tcga_rec
    
    
  } else{
    test_data <- test_data_survival
    analysis <- "survival"
    proof_genes_pt_tcga <- proof_genes_pt_tcga_survival
    time_roc <- time_roc_surv
    res_auc_tcga <- res_auc_tcga_surv
    
  }
  
  
  test_data <-
    test_data %>%
    mutate(EVENT_STAT = factor(EVENT_STAT))
  
  global_roc <- roc_curve(test_data, EVENT_STAT, risk_score) %>%
    mutate(label = "METABRIC")
  
  
  # 7.2 Loop that creates data frame with true positive and falsa positives of each time point
  
  plot_roc <- map_df(c(12, 36, 60, 72, 120), function(i) {
    data.frame(
      FP = time_roc$FP[, paste0("t=", i)],
      
      TP = time_roc$TP[, paste0("t=", i)],
      
      Time = factor(i),
      
      data_set = "METABRIC"
    )
    
  })
  
  
  # 7.3.2 Object with labels for the plot with the numerical values of the AUCs
  
  # 7.3.3 Same thing but for facet wrap labels
  
  facet_labels <-
    data.frame(
      Time = factor(c(12, 36, 60, 72, 120)),
      AUC_Text = paste0("AUC: ", round(100 * as.numeric(time_roc$AUC[1:5]), 3), "%"),
      data_set = "METABRIC"
    )
  
  
  
  # TCGA --------------------------------------------------------------------
  
  
  # 3.- ROC and AUC ---------------------------------------------------------
  
  # 3.1.2 View the AUC values
  
  res_auc_res_tcga <- res_auc_tcga$AUC %>%
    as.data.frame()
  
  # 3.1.3
  
  proof_genes_pt_tcga <-
    proof_genes_pt_tcga %>%
    mutate(EVENT_STAT = factor(EVENT_STAT))
  
  # 3.1.4 Global ROC curve object
  
  global_roc_tcga <- roc_curve(proof_genes_pt_tcga, EVENT_STAT, .pred_linear_pred) %>%
    mutate(label = "TCGA")
  
  
  # 3.2 Combine all time points into a long dataframe
  
  plot_roc_tcga <- map_df(c(12, 36, 60, 72, 120), function(i) {
    # This functions as a for loop
    time_label <- paste0("t=", i)
    
    data.frame(
      FP = res_auc_tcga$FP[, time_label],
      TP = res_auc_tcga$TP[, time_label],
      Time = factor(i),
      data_set = "TCGA"
    )
  })
  
  
  # 3.3 Labels for faceted plot
  
  facet_labels_tcga <- data.frame(
    Time = factor(c(12, 36, 60, 72, 120)),
    AUC_Text = paste0("AUC: ", 100 * round(as.numeric(
      res_auc_tcga$AUC[1:5]
    ), 3), "%"),
    data_set = "TCGA"
  )
  
  if (analysis == "recurrence") {
    
    global_roc_tcga_rec <- global_roc_tcga 
    plot_roc_tcga_rec <- plot_roc_tcga
    facet_labels_tcga_rec <- facet_labels_tcga
    
    global_roc_rec <- global_roc 
    plot_roc_rec <- plot_roc
    facet_labels_rec <- facet_labels
    
  } else{
    
    global_roc_tcga_surv <- global_roc_tcga
    plot_roc_tcga_surv <- plot_roc_tcga
    facet_labels_tcga_surv <- facet_labels_tcga
    
    global_roc_surv <- global_roc
    plot_roc_surv <- plot_roc
    facet_labels_surv <- facet_labels
    
  }
  
}


# gse2034 -----------------------------------------------------------------


# 2.3 View the AUC values
auc_gse2034 <- res_auc_gse2034$AUC
print(auc_gse2034)


proof_genes_pt_gse2034 <-
  proof_genes_pt_gse2034 %>%
  mutate(EVENT_STAT = factor(EVENT_STAT))

# 2.4 Object to later plot ROC curves

global_roc_gse2034 <- roc_curve(proof_genes_pt_gse2034, EVENT_STAT, .pred_linear_pred) %>%
  mutate(label = "GSE2034")


# 2.4.2 Combine all time points into a long dataframe

plot_roc_gse2034 <- map_df(c(12, 36, 60, 120), function(i) {
  # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc_gse2034$FP[, time_label],
    TP = res_auc_gse2034$TP[, time_label],
    Time = factor(i),
    data_set = "GSE2034"
  )
})


# 2.4.3 Labels for faceted plot

facet_labels_gse2034 <- data.frame(
  Time = factor(c(12, 36, 60, 120)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(
    res_auc_gse2034$AUC[1:4]
  ), 3), "%"),
  data_set = "GSE2034"
)


# GSE96058 ----------------------------------------------------------------




# 2.7 Combine all time points into a long dataframe

plot_roc_gse96058 <- map_df(c(12, 36, 60, 72), function(i) {
  # This functions as a for loop
  time_label <- paste0("t=", i)
  
  data.frame(
    FP = res_auc_gse96058$FP[, time_label],
    TP = res_auc_gse96058$TP[, time_label],
    Time = factor(i),
    data_set = "GSE96058"
  )
})


# 2.7.5 Global ROC object


proof_genes_pt_gse96058 <-
  proof_genes_pt_gse96058 %>%
  mutate(EVENT_STAT = factor(EVENT_STAT))

global_roc_gse96058 <- roc_curve(proof_genes_pt_gse96058, EVENT_STAT, .pred_linear_pred) %>%
  mutate(label = "GSE96058")

# 2.8 Extract AUC and SE for each time point

auc_table_gse96058 <- data.frame(
  Time = c(12, 36, 60, 72),
  `AUC (%)` = round(as.numeric(res_auc_gse96058$AUC[c("t=36", "t=60", "t=72", "t=80")]), 4),
  `SE` = res_auc_gse96058$times
)


# 2.8.2 Labels

legend_labels_gse96058 <- paste0("t=",
                                 auc_table_gse96058$Time,
                                 " (AUC: ",
                                 100 * auc_table_gse96058$AUC...,
                                 "%)")

# 2.8.3 Labels for faceted plot

facet_labels_gse96058 <- data.frame(
  Time = factor(c(12, 36, 60, 72)),
  AUC_Text = paste0("AUC: ", 100 * round(as.numeric(
    res_auc_gse96058$AUC[1:4]
  ), 3), "%"),
  data_set = "GSE96058"
)






# Object with all the ROC data from thje different databases

plot_roc_total_rec <- bind_rows(plot_roc_rec, plot_roc_gse2034, plot_roc_tcga_rec)

# Object with all the label data from all the databases

facet_labels_final_rec <- bind_rows(facet_labels_rec, facet_labels_tcga_rec, facet_labels_gse2034) %>%
  mutate(y_pos = case_when(
    data_set == "METABRIC" ~ 0.10,
    data_set == "GSE96058" ~ 0.18,
    data_set == "GSE2034" ~ 0.18,
    data_set == "TCGA"     ~ 0.26
  ))

# Plot

rec_1 <- ggplot(plot_roc_total_rec, aes(x = FP, y = TP, color = data_set)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap( ~ Time) +
  scale_color_manual(values = c("#D8AEDDFF", "#939aef", "#CB74ADFF")) +
  geom_text(data = facet_labels_final_rec, 
            aes(x = 0.75, y = y_pos, label = AUC_Text, color = data_set), 
            size = 5, fontface = "bold") +
  labs( 
    subtitle = "Recurrence",
    x = "1 - Specificity (FP)",
    y = "Sensitivity (TP)",
    color = "Cohort") +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 20, face = "bold"),  
    plot.subtitle = element_text(size = 16),              
    axis.title = element_text(size = 16),                 
    axis.text = element_text(size = 12),                  
    strip.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12),               
    legend.text = element_text(size = 12)                 
  )

########################################################################################


global_roc_all_rec <- bind_rows(global_roc_rec, global_roc_gse2034, global_roc_tcga_rec)


rec_glob <- ggplot(global_roc_all_rec, aes(x = 1 - specificity, y = sensitivity, color = label)) +
  geom_path(linewidth = 1.2) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_classic() +
  labs(title = "Global ROC curves for recurrence") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  scale_color_manual(values = c("#D8AEDDFF", "#939aef", "#CB74ADFF")) 








# Object with all the ROC data from thje different databases

plot_roc_total_surv <- bind_rows(plot_roc_surv, plot_roc_gse96058, plot_roc_tcga_surv)

# Object with all the label data from all the databases

facet_labels_final_surv <- bind_rows(facet_labels_surv, facet_labels_tcga_surv, facet_labels_gse96058) %>%
  mutate(y_pos = case_when(
    data_set == "METABRIC" ~ 0.10,
    data_set == "GSE96058" ~ 0.18,
    data_set == "GSE2034" ~ 0.18,
    data_set == "TCGA"     ~ 0.26
  ))

# Plot

surv_1 <- ggplot(plot_roc_total_surv, aes(x = FP, y = TP, color = data_set)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap( ~ Time) +
  scale_color_manual(values = c("#D8AEDDFF", "#939aef", "#CB74ADFF")) +
  geom_text(data = facet_labels_final_surv, 
            aes(x = 0.75, y = y_pos, label = AUC_Text, color = data_set), 
            size = 5, fontface = "bold") +
  labs( 
    title = "Time-Dependent ROC Curves",
    subtitle = "Survival analysis",
    x = "1 - Specificity (FP)",
    y = "Sensitivity (TP)",
    color = "Cohort") +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 20, face = "bold"),  
    plot.subtitle = element_text(size = 16),              
    axis.title = element_text(size = 16),                 
    axis.text = element_text(size = 12),                  
    strip.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12),               
    legend.text = element_text(size = 12)                 
  )

########################################################################################


global_roc_all_surv <- bind_rows(global_roc_surv, global_roc_gse96058, global_roc_tcga_surv)


surv_glob <- ggplot(global_roc_all_surv, aes(x = 1 - specificity, y = sensitivity, color = label)) +
  geom_path(linewidth = 1.2) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_classic() +
  labs(title = "Global ROC curves for survival",
       color = "Cohort") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold")
  ) + 
  scale_color_manual(values = c("#D8AEDDFF", "#939aef", "#CB74ADFF")) 
roc_plots <- 
  ((surv_glob | surv_1)/
    (rec_glob | rec_1)) 