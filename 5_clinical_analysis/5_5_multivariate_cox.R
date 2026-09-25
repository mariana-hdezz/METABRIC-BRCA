library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(flextable)
library(stringr)
library(broom)
library(patchwork)


# Load

test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
proof_genes_pt_cox_survival <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")

proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")
proof_genes_pt_gse2034_cox <- readRDS("./output_data/proof_genes_pt_gse2034_cox.RDS")

proof_genes_pt_gse96058  <- readRDS("./output_data/proof_genes_pt_gse96058.RDS")
proof_genes_pt_gse96058_cox  <- readRDS("./output_data/proof_genes_pt_gse96058_cox.RDS")

proof_genes_pt_tcga_cox_survival   <- readRDS("./output_data/proof_genes_pt_tcga_cox_survival.RDS")
proof_genes_pt_tcga_cox_recurrence <- readRDS("./output_data/proof_genes_pt_tcga_cox_recurrence.RDS")

for (i in c(1:2)) {
  if (i == 1) {
    analysis <- "recurrence"
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_recurrence
    proof_genes_pt_cox <- proof_genes_pt_cox_recurrence
    
    
  } else{
    analysis <- "survival"
    proof_genes_pt_tcga_cox <- proof_genes_pt_tcga_cox_survival
    proof_genes_pt_cox <- proof_genes_pt_cox_survival
    
  }
  # 8.- Independence test ---------------------------------------------------
  
  # 8.3 Multivariate cox comparing the variables including the score to test its independent value
  
  cox_model <- coxph(
    surv_obj ~ NPI + HORMONE + CHEMO + RADIO + strata(SURGERY) + MENO + strata(HER2) + AGE + SCORE +  PAM50 + strata(INTCLUST) + strata(HIST),
    data = proof_genes_pt_cox
  )
  
  summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, data = proof_genes_pt_cox))
  
  
  summary(cox_model)
  
  # 8.3.2 Tidy format
  
  independent_prog <- cox_model %>%
    tidy(exponentiate = TRUE, conf.int = TRUE)
  
  supplementary_table <- independent_prog %>%
    mutate(
      Cohort = "METABRIC",
      Feature = recode_values(
        term,
        "SCORE"              ~ "Signature Score",
        "LYMPH"              ~ "Lymph Node Status",
        "AGE"                ~ "Age at Diagnosis",
        "NPI"                ~ "Nottingham Prognostic Index",
        "CHEMOYES"           ~ "Chemotherapy (Yes)",
        "HORMONEYES"         ~ "Hormone Therapy (Yes)",
        "MENOPre"            ~ "Menopausal Status (Pre)",
        "PAM50LumA"          ~ "PAM50 Luminal A",
        "PAM50LumB"          ~ "PAM50 Luminal B",
        "PAM50Her2"          ~ "PAM50 Her2-enriched",
        "PAM50claudin-low"   ~ "PAM50 Claudin-low",
        "PAM50Normal"        ~ "PAM50 Normal-like",
        "PAM50NC"        ~ "PAM50 Not classified",
        "RADIOYES" ~ "Radiotherapy (Yes)",
        default = term
      ),
      
      `Hazard Ratio (HR)` = round(estimate, 3),
      `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
      
      `p-value` = scales::scientific(p.value, digits = 3)
    ) %>%
    dplyr::select(Cohort, Feature,
                  `Hazard Ratio (HR)`,
                  `95% Confidence Interval`,
                  `p-value`) %>%
    arrange(`Hazard Ratio (HR)`)
  
  
  # 9.- Results -------------------------------------------------------------
  
  # 9.1 Index number of the parameters to print the evaluation and to graph
  
  num_param_compare <- c(1, 2, 3, 4, 6, 7, 10)
  
  
  # 9.3 Forest plot of the evaluated parameters
  
  
  cox_p_metabric <- independent_prog[num_param_compare, ] %>%
    filter(estimate > 0.0001, conf.high < 100) %>%
    mutate(
      term = recode_values(
        term,
        "SCORE"              ~ "Signature Score",
        "LYMPH"              ~ "Lymph Node Status",
        "AGE"                ~ "Age at Diagnosis",
        "NPI"                ~ "Nottingham Prognostic Index",
        "CHEMOYES"           ~ "Chemotherapy (Yes)",
        "HORMONEYES"         ~ "Hormone Therapy (Yes)",
        "SURGERYMASTECTOMY"  ~ "Mastectomy",
        "PAM50LumA"          ~ "Claudin subtype Luminal A",
        "PAM50LumB"          ~ "Claudin subtype Luminal B",
        "PAM50Her2"          ~ "Claudin subtype Her2-enriched",
        "PAM50claudin-low"   ~ "Claudin subtype Claudin-low",
        "PAM50Normal"        ~ "Claudin subtype Normal-like",
        "HISTMucinous"       ~ "Histology: Mucinous",
        "HISTMixed"          ~ "Histology: Mixed",
        "HISTLobular"        ~ "Histology: Lobular",
        "HISTMedullary"      ~ "Histology: Medullary",
        "INTCLUST2"          ~ "IntClust 2",
        "INTCLUST3"          ~ "IntClust 3",
        "INTCLUST5"          ~ "IntClust 5",
        "INTCLUST6"          ~ "IntClust 6",
        "INTCLUST7"          ~ "IntClust 7",
        "INTCLUST8"          ~ "IntClust 8",
        "INTCLUST9"          ~ "IntClust 9",
        "MENOPre" ~ "Premenopause",
        "RADIOYES" ~ "Radiotherapy (Yes)",
        default = term
      ),
      term = reorder(term, estimate),
      significant = p.value < 0.05
    ) %>%
    ggplot(aes(x = estimate, y = term, color = significant)) +
    geom_point() +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  width = 0.5,
                  linewidth = 1.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    labs(x = "Hazard Ratio (log scale)", y = "Clinical & Molecular Features", color = "Significance (p < 0.05)", title = paste0("Multivariate Cox: ", str_to_title(analysis), " METABRIC")) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",) +
    scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) +
    labs(
      x = "Hazard Ratio (log scale)",
      y = "",
      color = "Significance (p < 0.05)",
    )
  
  
  
  # TCGA --------------------------------------------------------------------
  
  
  # 4.- Multivariate cox
  
  
  # 4.2 Actual cox model with parameters to evaluare
  
  cox_model_tcga <- coxph(
    surv_obj ~ SCORE + AGE + LYMPH + strata(PAM50) + strata(HER2) + strata(HIST) + RADIO + strata(NEO) + TARG_TX,
    data = proof_genes_pt_tcga_cox
  )
  
  # 4.3 Tidy format
  
  independent_prog_tcga <- cox_model_tcga %>%
    tidy(exponentiate = TRUE, conf.int = TRUE)
  
  supplementary_table_tcga <- independent_prog_tcga %>%
    mutate(
      Cohort = "TCGA",
      Feature = recode_values(
        term,
        "SCORE"              ~ "Signature Score",
        "LYMPH"              ~ "Lymph Node Status",
        "AGE"                ~ "Age at Diagnosis",
        "RADIOYES" ~ "Radiotherapy",
        "NEOYes" ~ "Neoadjuvant chemotherapy",
        "TARG_TXYES" ~ "Targeted treatment",
        default = term
      ),
      
      `Hazard Ratio (HR)` = round(estimate, 3),
      `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
      
      `p-value` = scales::scientific(p.value, digits = 3)
    ) %>%
    dplyr::select(Cohort, Feature,
                  `Hazard Ratio (HR)`,
                  `95% Confidence Interval`,
                  `p-value`) %>%
    arrange(`Hazard Ratio (HR)`)
  
  flextable(supplementary_table_tcga) %>%
    autofit()
  
  # 4.3.2 Results
  
  summary(cox_model_tcga)
  
  # 4.4 Results to paragraph
  
  num_param_compare <- c(1:9)
  
  # 4.5 Forest plot ignoring values that tend to infinite
  
  cox_p_tcga <- independent_prog_tcga %>%
    filter(estimate > 0.0001, conf.high < 100) %>%
    mutate(
      term = recode_values(
        term,
        "SCORE"      ~ "Signature Score",
        "LYMPH"      ~ "Lymph Node Status",
        "AGE"        ~ "Age at Diagnosis",
        "RADIOYES"   ~ "Radiotherapy",
        "NEOYes"     ~ "Neoadjuvant chemotherapy",
        "TARG_TXYES" ~ "Targeted treatment",
        default = term
      ),
      term = reorder(term, p.value),
      significant = factor(p.value < 0.05, levels = c("FALSE", "TRUE"))
    ) %>%
    ggplot(aes(x = estimate, y = term, color = significant)) +
    
    # FIX: Add show.legend = TRUE to both geoms to force ggplot to draw the keys
    geom_point(show.legend = TRUE) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.5,
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(
        color = "black",
        fill = "white",
        linewidth = 0.5
      ),
      legend.box.background = element_rect(color = "black", linewidth = 1),
      legend.key = element_rect(color = "gray80", linewidth = 0.5),
      legend.position = ifelse(analysis == "recurrence","right" , "none")
    ) +
    scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD"),
                       drop = FALSE) +
    labs(
      x = "Hazard Ratio (log scale)",
      y = "Clinical & Molecular Features",
      color = "Significance (p < 0.05)",
      title = paste0("Multivariate Cox: ", str_to_title(analysis), " TCGA")
    )
  
  
  if (analysis == "recurrence") {
    supplementary_table_tcga_rec <- supplementary_table_tcga
    cox_p_tcga_rec <- cox_p_tcga 
    supplementary_table_rec <- supplementary_table
    cox_p_metabric_rec <- cox_p_metabric
    saveRDS(cox_model_tcga, "./output_data/cox_model_tcga_rec.RDS")
    saveRDS(cox_model, "./output_data/cox_model_rec.RDS")
      
  } else{
    supplementary_table_tcga_surv <- supplementary_table_tcga
    cox_p_tcga_surv <- cox_p_tcga
    supplementary_table_surv <- supplementary_table
    cox_p_metabric_surv <- cox_p_metabric
    saveRDS(cox_model_tcga, "./output_data/cox_model_tcga_surv.RDS")
    saveRDS(cox_model, "./output_data/cox_model_surv.RDS")
  }
  
}

# GSE2034 -----------------------------------------------------------------

# 3.1.2 Once again split to test the cox with score as a continuous variable and not as a divided low and high risk

gse2034_split_cox <- survSplit(
  formula = Surv(EVENT_MON, EVENT_STAT) ~ .,
  data = proof_genes_pt_gse2034_cox,
  cut = 65,
  episode = "time_group",
  id = "patient_id"
)

# 3.2 Fit the model with the interaction

cox_model_gse2034 <- coxph(Surv(tstart, EVENT_MON, EVENT_STAT) ~ SCORE:strata(time_group),
                           data = gse2034_split_cox)

# 3.3 Tidy

independent_prog_gse2034 <-
  cox_model_gse2034 %>%
  tidy(exponentiate = TRUE, conf.int = TRUE)

# 3.3.2 Clean names for table

supplementary_table_gse2034 <- independent_prog_gse2034 %>%
  mutate(
    Cohort = "GSE2034",
    Feature = recode_values(
      term,
      "SCORE:strata(time_group)time_group=1"              ~ "Time Group 1 (>65 months)",
      "SCORE:strata(time_group)time_group=2"              ~ "Time Group 2 (<65 months)",
      default = term
    ),
    
    `Hazard Ratio (HR)` = round(estimate, 3),
    `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
    
    `p-value` = scales::scientific(p.value, digits = 3)
  ) %>%
  dplyr::select(Cohort, Feature,
                `Hazard Ratio (HR)`,
                `95% Confidence Interval`,
                `p-value`) %>%
  arrange(`Hazard Ratio (HR)`)

# 3.3.3 Table

flextable(supplementary_table_gse2034) %>%
  autofit()


# 3.4 Forest plot ignoring values that tend to infinite

cox_p_gse2034 <-  independent_prog_gse2034 %>%
  filter(estimate > 0.0001, conf.high < 100) %>%
  mutate(
    term = recode_values(
      term,
      "SCORE:strata(time_group)time_group=1"              ~ "Time Group 1 (>65 months)",
      "SCORE:strata(time_group)time_group=2"              ~ "Time Group 2 (<65 months)",
      default = term
    ),
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbar(
    orientation = "y",
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.5,
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_linedraw() +
  theme_classic(base_size = 15) +
  ggtitle("Multivariate Cox: Recurrence GSE2034") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD")) +
  labs(x = "Hazard Ratio (log scale)",
       y = "",
       color = "Significance (p < 0.05)")



# GSE96058 ----------------------------------------------------------------


# 3.- Independence test ---------------------------------------------------


summary(coxph(surv_obj ~ AGE + LYMPH + SCORE, data = proof_genes_pt_gse96058_cox))

# 3.2 Multivariate cox

cox_model_gse96058 <- coxph(surv_obj ~ PAM50 + KI67 + HER2 + AGE + LYMPH + CHEMO + HORMONE + SCORE,
                            data = proof_genes_pt_gse96058_cox)

independent_prog_gse96058 <- cox_model_gse96058 %>%
  tidy(exponentiate = TRUE, conf.int = TRUE)

summary(cox_model_gse96058)


# 4.- Results -------------------------------------------------------------

# 4.1 Index numbers to use of the multivariate cox

num_param_compare_gse <- c(2, 3, 5, 7, 8, 10, 11, 12, 13, 14)


supplementary_table_gse96058 <- independent_prog_gse96058 %>%
 mutate(
   Cohort = "GSE96058",
    Feature = recode_values(
      term,
      "SCORE"              ~ "Signature Score",
      "LYMPH4toX"              ~ "Lymph nodes >= 4",
      "LYMPHNodeNegative"              ~ "Lymph nodes negative",
      "LYMPHSubMicroMet"              ~ "Lymph micrometastasis",
      "AGE"                ~ "Age at Diagnosis",
      "CHEMO"           ~ "Chemotherapy",
      "HORMONE"         ~ "Hormone Therapy",
      "PAM50LumA"          ~ "Claudin subtype Luminal A",
      "PAM50LumB"          ~ "Claudin subtype Luminal B",
      "PAM50Her2"          ~ "Claudin subtype Her2-enriched",
      "PAM50Normal"        ~ "Claudin subtype Normal-like",
      default = term
    ),
    
    `Hazard Ratio (HR)` = round(estimate, 3),
    `95% Confidence Interval` = paste0(round(conf.low, 3), " – ", round(conf.high, 3)),
    
    `p-value` = scales::scientific(p.value, digits = 3)
  ) %>%
  dplyr::select(Cohort, Feature,
                `Hazard Ratio (HR)`,
                `95% Confidence Interval`,
                `p-value`) %>%
  arrange(`Hazard Ratio (HR)`)

# 4.3 Forest Plot

cox_p_gse96058 <- independent_prog_gse96058[num_param_compare_gse, ] %>%
  filter(estimate > 0.0001, conf.high < 100) %>%
  filter(!(term == "LYMPHNA")) %>%
  mutate(
    term = recode_values(
      term,
      "SCORE"              ~ "Signature Score",
      "LYMPH4toX"              ~ "Lymph nodes >= 4",
      "LYMPHNodeNegative"              ~ "Lymph nodes negative",
      "LYMPHSubMicroMet"              ~ "Lymph micrometastasis",
      "AGE"                ~ "Age at Diagnosis",
      "CHEMO"           ~ "Chemotherapy",
      "HORMONE"         ~ "Hormone Therapy",
      "PAM50LumA"          ~ "Claudin subtype Luminal A",
      "PAM50LumB"          ~ "Claudin subtype Luminal B",
      "PAM50Her2"          ~ "Claudin subtype Her2-enriched",
      "PAM50Normal"        ~ "Claudin subtype Normal-like",
      default = term
    ),
    term = reorder(term, estimate),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.5,
                 linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)",
       y = "",
       color = "Significance (p < 0.05)") +
  theme_classic(base_size = 15) +
  ggtitle("Multivariate Cox: Survival GSE96058") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
  ) +
  scale_color_manual(values = c("FALSE" = "#68228b", "TRUE" = "#3477FD"))


saveRDS(cox_model_gse96058, "./output_data/cox_model_gse96058.RDS")
saveRDS(cox_model_gse2034, "./output_data/cox_model_gse2034.RDS")

((cox_p_metabric_surv / cox_p_tcga_surv / cox_p_gse96058) |
(cox_p_metabric_rec / cox_p_tcga_rec / cox_p_gse2034) ) + 
  plot_annotation(
    tag_levels = "A"
  )




rbind(
"Survival",
supplementary_table_surv,
supplementary_table_gse96058,
supplementary_table_tcga_surv
) %>% 
  flextable() %>% 
  merge_at(i = 1, j = 1:5) %>% 
  merge_at(i = 2:13, j = 1) %>% 
  merge_at(i = 14:27, j = 1) %>% 
  merge_at(i = 28:32, j = 1) %>% 
  hline(i = 13) %>% 
  hline(i = 27) 


rbind(
"Recurrence",
supplementary_table_rec,
supplementary_table_tcga_rec,
supplementary_table_gse2034
)  %>% 
  flextable() %>% 
  merge_at(i = 1, j = 1:5) %>% 
  merge_at(i = 2:14, j = 1) %>% 
  merge_at(i = 15:19, j = 1) %>% 
  merge_at(i = 20:21, j = 1)  %>% 
  hline(i = 14) %>% 
  hline(i = 19) 
