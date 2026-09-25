library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(purrr)
library(flextable)

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
  
  
  # 6.4.0 Convert EVENT_sTAT to factor
  
  proof_genes_pt_cox <-
    proof_genes_pt_cox %>%
    mutate(EVENT_STAT = factor(EVENT_STAT))
  
  # 6.4.3 Prepare to plot the different comparisons
  
  proof_genes_pt_long <-
    proof_genes_pt_cox %>%
    pivot_longer(
      cols = c(SURGERY, CHEMO, HORMONE),
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    mutate(
      Value = case_when(
        toupper(Value) == "YES" ~ "Yes",
        toupper(Value) == "NO"  ~ "No",
        Value == "BREAST CONSERVING" ~ "Breast conserving",
        Value == "MASTECTOMY" ~ "Matectomy",
        TRUE ~ Value
      )
    )
  
  
  # 6.5 Wilcox test between treatment types
  
  wilcox_treatment <-
    proof_genes_pt_long %>%
    dplyr::select(EVENT_STAT, Parameter, Value, SCORE) %>%
    group_by(Parameter, Value) %>%
    filter(n() != 1) %>%
    group_modify(~ {
      test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
      tidy(test)
    }) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p.value, method = "holm"), )
  
  # 6.5.2 Cox analysis of interaction between score and treatment
  
  tx_vars <- c("CHEMO", "HORMONE", "SURGERY")
  
  cox_sum <- list()
  
  for (t in tx_vars) {
    formula <- as.formula(paste("surv_obj ~ ", t, " * SCORE")) # Establish formula
    
    proof_genes_pt_txcox <-
      proof_genes_pt_cox %>%
      group_by(.data[[t]]) %>%
      filter(
        n_distinct(EVENT_STAT) == 2,
        # FIlter data with complete separation!(is.na(.data[[t]])),
        # Filter NA!(.data[[t]] == "") # Filter fod the unregistered surgery
      ) %>%
      ungroup()
    
    cox_sum[[paste(i, t)]] <- summary(coxph(formula , data = proof_genes_pt_txcox))
    print(paste(i, t))
    
  }
  
  
  # 6.5.3 Plot
  
  score_tx <- ggplot(data = proof_genes_pt_long, aes(y = SCORE, x = Value, fill = Value)) +
    geom_boxplot() +
    facet_wrap(
      ~ Parameter + EVENT_STAT,
      scales = "free_x",
      ncol = 2,
      labeller = labeller(
        Parameter = c(
          "CHEMO" = "Chemotherapy",
          "HORMONE" = "Hormone therapy",
          "SURGERY" = "Surgery Modality"
        ),
        EVENT_STAT = c(
          "0" = ifelse(analysis == "recurrence" , "Non-recurrent", "Alive"),
          "1" = ifelse(analysis == "recurrence" , "Recurrent", "Deceased")
        )
      )
    ) +
    scale_fill_paletteer_d("khroma::iridescent", direction = -1) +
    labs(title = paste0(
      "Score change on treatment: METABRIC (",
      stringr::str_to_title(analysis),
      ")"
    )) +
    theme_classic(base_size = 12) +
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
    labs(x = "")
  
  
  
  # TCGA --------------------------------------------------------------------
  
  # treatment comparisons
  
  vars <- c("RADIO", "NEO", "OTHER_TX", "TARG_TX")
  
  # To make data frame do the "loop" with map_dfr
  
  wilcox_treatment_tcga <- map_dfr(vars, function(var) {
    df <- proof_genes_pt_tcga_cox %>%
      group_by(.data[[var]]) %>%
      filter(n_distinct(EVENT_STAT) == 2) %>%
      group_modify(~ {
        test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
        tidy(test)
      }) %>%
      mutate(level = cur_group()[[1]])  %>%
      ungroup() %>%
      mutate(variable = var,
             adj_p_value = p.adjust(p.value, method = "holm")) %>%
      dplyr::select(variable, level, p.value, adj_p_value)
  })
  
  
  
  # 5.3 Pivot longer so as to plot the different treatment given scores faceted by event status
  
  proof_genes_pt_tcga_long <-
    proof_genes_pt_tcga_cox %>%
    pivot_longer(
      cols = c(NEO, SURGERY, RADIO, TARG_TX),
      # Add any other parameters here
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    mutate(
      Value = case_when(
        toupper(Value) == "YES" ~ "Yes",
        toupper(Value) == "NO"  ~ "No",
        TRUE ~ Value # Keeps things like "Lumpectomy" and "NA" untouched
      ),
      Value = case_when(
        tolower(Value) == "mastectomy nos" ~ "Mastectomy",
        tolower(Value) == "modified radical mastectomy" ~ "Modified Radical Mastectomy",
        TRUE ~ stringr::str_to_title(Value) # Capitalizes first letter of other categories like "Lumpectomy", "Other"
      )
    )
  
  # 5.3.2 Plot
  
  
  score_tx_tcga <- ggplot(data = proof_genes_pt_tcga_long %>% drop_na(Value), aes(y = SCORE, x = Value, fill = Value)) +
    geom_boxplot() +
    facet_wrap(
      ~ Parameter + EVENT_STAT,
      scales = "free_x",
      ncol = 2,
      labeller = labeller(
        EVENT_STAT = c(
          "0" = ifelse(analysis == "recurrence" , "Non-recurrent", "Alive"),
          "1" = ifelse(analysis == "recurrence" , "Recurrent", "Deceased")
        ),
        Parameter  = as_labeller(
          c(
            "NEO" = "Neoadjuvant therapy",
            "RADIO" = "Radiotherapy",
            "SURGERY" = "Surgery modality",
            "TARG_TX" = "Targeted treatment"
          )
        )
      )
    ) +
    scale_fill_paletteer_d(
      "khroma::iridescent",
      direction = -1,
      labels = function(x)
        stringr::str_to_title(x)
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = "Treatment modality",
      y = "",
      fill = "Treated",
      title = paste0(
        "Score change on treatment: TCGA (",
        stringr::str_to_title(analysis),
        ")"
      )
    ) +
    scale_x_discrete(labels = c("0" = "Untreated", "1" = "Treated")) +
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
      axis.text.x = element_text(
        angle = 22,
        hjust = 1,
        vjust = 1
      )
    )
  
  # 5.4 Cox based on treatment
  vars_strat <- c("RADIO", "OTHER_TX", "TARG_TX")
  
  cox_sum_tcga <- list()
  
  for (t in seq_along(vars_strat)) {
    formula <- as.formula(paste("surv_obj ~ ", vars_strat[t], " * SCORE")) # Establish formula
    
    proof_genes_pt_tcga_wilcox <-
      proof_genes_pt_tcga_cox %>%
      group_by(.data[[vars_strat[t]]]) %>%
      filter(n_distinct(EVENT_STAT) == 2, ) %>%
      ungroup()
    
    cox_sum_tcga[[paste(i, vars_strat[t])]] <- summary(coxph(formula , data = proof_genes_pt_tcga_wilcox))
    print(paste(i, t))
    
  }
  
  print(i)
  
  if (analysis == "recurrence") {
    wilcox_treatment_tcga_rec <- wilcox_treatment_tcga
    wilcox_treatment_rec <- wilcox_treatment
    score_tx_tcga_rec <- score_tx_tcga
    score_tx_rec <- score_tx
    tx_cox_txga_rec <- cox_sum_tcga
    tx_cox_rec <- cox_sum
    
  } else{
    wilcox_treatment_tcga_surv <- wilcox_treatment_tcga
    wilcox_treatment_surv <- wilcox_treatment
    score_tx_tcga_surv <- score_tx_tcga
    score_tx_surv      <- score_tx
    tx_cox_txga_surv <- cox_sum_tcga
    tx_cox_surv <- cox_sum
    
  }
  
}




wilcox_treatment_gse96058 <- list()

# 4.5 Wilcoxon test for treatment

for (t in tx_vars[1:2]) {
  wilcox_treatment_gse96058[[t]] <-
    proof_genes_pt_gse96058_cox %>%
    dplyr::select(EVENT_STAT, all_of(t), SCORE) %>%
    group_by(.data[[t]]) %>%
    group_modify( ~ {
      test <- wilcox.test(SCORE ~ EVENT_STAT, data = .)
      tidy(test)
    }) %>%
    ungroup() %>%
    mutate(adj_p_value = p.adjust(p.value, method = "holm"))
}

# 4.6.2 Long format to visualize how the scorescores patient based on event an treatment

proof_genes_long_gse96058 <- proof_genes_pt_gse96058_cox %>%
  mutate(HORMONE = factor(HORMONE), CHEMO = factor(CHEMO)) %>%
  pivot_longer(
    cols = c(HORMONE, CHEMO),
    names_to = "Parameter",
    values_to = "Value"
  )


# 4.6.3  Plot

score_tx_gse96058 <- ggplot(proof_genes_long_gse96058,
                            aes(
                              y = SCORE,
                              x = Value,
                              fill = factor(Value, labels = c("0" = "Untreated", "1" = "Treated"))
                            )) +
  geom_boxplot() +
  facet_wrap(
    ~ Parameter + EVENT_STAT,
    scales = "free_x",
    labeller = labeller(EVENT_STAT = as_labeller(c(
      "0" = "Alive", "1" = "Deceased"
    )), Parameter  = as_labeller(
      c("CHEMO" = "Chemotherapy", "HORMONE" = "Hormonal treatment")
    ))
  ) +
  scale_fill_paletteer_d("khroma::iridescent", direction = -1) +
  theme_classic(base_size = 12) +
  ggtitle("Score change on treatment: GSE96058 (Survival)") +
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
  labs(y = "", x  = "")

tx <- c("HORMONE", "CHEMO")

tx_cox_gse96058 <- list()

for (i in seq_along(tx)) {
  formula_gse96058 <- as.formula(paste("surv_obj ~ ", tx[i], " * SCORE", collapse = "")) # Establish formula
  
  tx_cox_gse96058[[paste(i, tx[i])]] <- summary(coxph(formula_gse96058, data = proof_genes_pt_gse96058_cox))
  
}




(((score_tx_surv |
     score_tx_tcga_surv |
     score_tx_gse96058) |
    (score_tx_rec |
       score_tx_tcga_rec)
) +
    plot_layout(widths = c(1, 1, 1, 2))) +
  plot_annotation(tag_levels = "A")


colapse_list <- lapply(seq_along(wilcox_treatment_gse96058), function(i) {
  wilcox_treatment_gse96058[[i]] %>%
    as.data.frame() %>%
    pivot_longer(cols = 1,
                 names_to = "Variable",
                 values_to = "Value")
})


wilcox_gse96058_df <- do.call(rbind, colapse_list)





tx_cox_surv
tx_cox_txga_surv
tx_cox_gse96058

rbind(
  "Survival",
  wilcox_treatment_surv %>%
    dplyr::select(Parameter, Value, adj_p_value) %>%
    rename("Variable" = Parameter, "Adjusted p value" = adj_p_value) %>%
    mutate(Cohort = "METABRIC")
  ,
  
  wilcox_treatment_tcga_surv %>%
    dplyr::select(variable, level, adj_p_value) %>%
    rename(
      "Value" = level,
      "Adjusted p value" = adj_p_value,
      "Variable" = variable
    ) %>%
    mutate(Cohort = "TCGA"),
  
  wilcox_gse96058_df %>%
    dplyr::select(Variable, Value, adj_p_value) %>%
    rename("Adjusted p value" = adj_p_value) %>%
    mutate(Cohort = "GSE96058"),
  
  "Recurrence",
  wilcox_treatment_rec %>%
    dplyr::select(Parameter, Value, adj_p_value) %>%
    rename("Variable" = Parameter, "Adjusted p value" = adj_p_value) %>%
    mutate(Cohort = "METABRIC"),
  
  wilcox_treatment_tcga_rec %>%
    dplyr::select(variable, level, adj_p_value) %>%
    rename(
      "Value" = level,
      "Adjusted p value" = adj_p_value,
      "Variable" = variable
    ) %>%
    mutate(Cohort = "TCGA")
  
) %>%
  filter(! Value == "") %>% 
  relocate(Cohort) %>% 
  mutate(Variable = case_when(
    Variable == "CHEMO" ~ "Chemotherapy",
    Variable == "HORMONE" ~ "Hormone therapy",
    Variable == "SURGERY" ~ "Surgical intervention",
    Variable == "RADIO" ~ "Radiotherapy",
    Variable == "NEO" ~ "Neoadjuvant chemotherapy",
    Variable == "OTHER_TX" ~ "Other treatment",
    Variable == "TARG_TX" ~ "Targeted treatment",
  )) %>% 
flextable() %>%
  merge_at(i = 1, j = 1:4) %>%
  merge_at(i = 20, j = 1:4) %>% 
  merge_at(i = 2:7, j = 1) %>% 
  hline(i = 7) %>% 
  merge_at(i = 8:15, j = 1) %>% 
  hline(i = 15) %>% 
  merge_at(i = 16:19, j = 1) %>% 
  hline(i = 19) %>% 
  merge_at(i = 21:26, j = 1) %>% 
  hline(i = 26) %>% 
  merge_at(i = 27:33, j = 1) 
  
  
score_tx_tcga_rec
score_tx_rec
tx_cox_txga_rec
tx_cox_rec
