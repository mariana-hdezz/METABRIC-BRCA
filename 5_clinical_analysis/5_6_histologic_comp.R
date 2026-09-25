library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(flextable)
library(paletteer)

# Load

test_data_survival <- readRDS("./output_data/test_data_survival.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
proof_genes_pt_cox_survival <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")

proof_genes_pt_gse2034 <- readRDS("./output_data/proof_genes_pt_gse2034.RDS")

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
  
  # METABRIC ----------------------------------------------------------------
  
  
  # 6.4 Wilcox test of Histology, Claudin subtype and INTCLUST data
  
  ch_var <- c("HIST", "INTCLUST", "PAM50")
  
  wilcoxon_comp <- list()
  
  for (t in ch_var) {
    wilcoxon_comp[[t]] <-
      proof_genes_pt_cox %>%
      dplyr::select(EVENT_STAT, all_of(t), SCORE) %>%
      filter(!is.na(.data[[t]])) %>%
      group_by(.data[[t]]) %>%
      filter(n_distinct(EVENT_STAT) == 2) %>%
      group_modify( ~ {
        test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
        tidy(test)
      }) %>%
      ungroup() %>%
      mutate(adj_p_value = p.adjust(p.value, method = "holm"))
    
  }
  
  # 6.4.2 Plot the  comparisons betweeen PAM50  subtypes
  
  score_subtype <- ggplot(proof_genes_pt_cox, aes(x = PAM50, y = SCORE, fill = PAM50)) +
    geom_boxplot() +
    scale_fill_paletteer_d("Redmonder::dPBIPuOr") +
    facet_wrap(~ EVENT_STAT, labeller = labeller(EVENT_STAT = c(
      "0" = ifelse(analysis == "recurrence", "Non-recurrent", "Alive"),
      "1" = ifelse(analysis == "recurrence", "Recurrent", "Deceased")
    ))) +
    theme_classic(base_size = 22) +
    ggtitle("Score change on subtype: METABRIC") +
    theme(
      strip.background = element_rect(
        fill = "black",
        color = "black",
        linewidth = 1
      ),
      strip.text = element_text(
        color = "white",
        face = "bold",
        size = 12
      ),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  # TCGA --------------------------------------------------------------------
  
  # 5.2 Wilcoxon test where on each PAM50 group we compare deceased vs live patients
  
  
  wilcoxon_comp_tcga <- list()
  
  for (t in ch_var) {
    wilcoxon_comp_tcga[[t]] <-
      proof_genes_pt_tcga_cox %>%
      dplyr::select(EVENT_STAT, all_of(t), SCORE) %>%
      filter(!is.na(.data[[t]])) %>%
      group_by(.data[[t]]) %>%
      filter(n_distinct(EVENT_STAT) == 2) %>%
      group_modify( ~ {
        test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
        tidy(test)
      }) %>%
      ungroup() %>%
      mutate(adj_p_value = p.adjust(p.value, method = "holm"))
  }
  
  # 5.1.2 Boxplot comparing to PAM50 facet wrapped by event status
  
  
  score_subtype_tcga <-
    ggplot(proof_genes_pt_tcga_cox %>% drop_na(PAM50),
           aes(y = SCORE, x = PAM50, fill = PAM50)) +
    geom_boxplot() +
    facet_wrap(~ EVENT_STAT,
               scales = "free_x",
               labeller = labeller(EVENT_STAT = c(
                 "0" = ifelse(analysis == "recurrence", "Non-recurrent", "Alive"),
                 "1" = ifelse(analysis == "recurrence", "Recurrent", "Deceased")
               ))) +
    theme_classic(base_size = 22) +
    ggtitle("Score change on subtype: TCGA") +
    theme(
      strip.background = element_rect(
        fill = "black",
        color = "black",
        linewidth = 1
      ),
      strip.text = element_text(
        color = "white",
        face = "bold",
        size = 12
      ),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    labs(y = "Score", x = "Intrinsic subtype") +
    scale_fill_paletteer_d("Redmonder::dPBIPuOr")
  
  
  if (analysis == "recurrence") {
    wilcoxon_comp_tcga_rec <- wilcoxon_comp_tcga
    wilcoxon_comp_rec <- wilcoxon_comp
    score_subtype_tcga_rec <- score_subtype_tcga
    score_subtype_rec <- score_subtype
    
    
  } else{
    wilcoxon_comp_tcga_surv <- wilcoxon_comp_tcga
    wilcoxon_comp_surv <- wilcoxon_comp
    score_subtype_tcga_surv <- score_subtype_tcga
    score_subtype_surv <- score_subtype
    
  }
  
}



# 4.5 Wilcoxon test for PAM50


wilcoxon_comp_gse96058 <-
  proof_genes_pt_gse96058_cox %>%
  dplyr::select(EVENT_STAT, PAM50, SCORE) %>%
  group_by(PAM50) %>%
  group_modify( ~ {
    test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
    tidy(test)
  }) %>%
  ungroup() %>%
  mutate(adj_p_value = p.adjust(p.value, method = "holm"))

# 4.4 Boxplot comparing to PAM50


score_subtype_gse96058 <- ggplot(proof_genes_pt_gse96058_cox,
                                 aes(y = SCORE, x = PAM50, fill = PAM50)) +
  geom_boxplot() +
  facet_wrap( ~ EVENT_STAT, labeller = labeller(EVENT_STAT = c("0" = "Alive", "1" = "Deceased"))) +
  theme_classic(base_size = 22) +
  ggtitle("Score change on subtype: GSE96058") +
  theme(
    strip.background = element_rect(
      fill = "black",
      color = "black",
      linewidth = 1
    ),
    strip.text = element_text(
      color = "white",
      face = "bold",
      size = 12
    ),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  labs(y = "", x = "Intrinsic Subtype") +
  scale_fill_paletteer_d("Redmonder::dPBIPuOr")

wilcoxon_comp_rec
wilcoxon_comp_tcga_rec 



wilcoxon_comp_surv 
wilcoxon_comp_tcga_surv
wilcoxon_comp_gse96058



(
(
  score_subtype_surv /
    score_subtype_tcga_surv /
    score_subtype_gse96058
)
 |
    (
      score_subtype_rec /
        score_subtype_tcga_rec
    )) +
    plot_annotation(tag_levels = "A")
