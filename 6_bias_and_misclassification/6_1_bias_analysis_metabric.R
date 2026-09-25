library(ggrepel)
library(paletteer)
library(patchwork)
library(tidymodels)
library(survminer)
library(dplyr)
library(tidyr)
library(tibble)

# load

final_fit_survival    <- readRDS("./output_data/final_fit_survival.RDS")
train_data_recurrence <- readRDS("./output_data/train_data_recurrence.RDS")
test_data_recurrence <- readRDS("./output_data/test_data_recurrence.RDS")
cox_model_recurrence  <- readRDS("./output_data/cox_model_rec.RDS")

final_fit_recurrence  <- readRDS("./output_data/final_fit_recurrence.RDS")
train_data_survival   <- readRDS("./output_data/train_data_survival.RDS")
test_data_survival   <- readRDS("./output_data/test_data_survival.RDS")
cox_model_survival    <- readRDS("./output_data/cox_model_surv.RDS")

proof_genes_pt_cox_survival <- readRDS("./output_data/proof_genes_pt_cox_survival.RDS")
proof_genes_pt_cox_recurrence <- readRDS("./output_data/proof_genes_pt_cox_recurrence.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    analysis <- "recurrence"
    final_fit  <- final_fit_recurrence
    train_data <- train_data_recurrence
    cox_model  <- cox_model_recurrence
    proof_genes_pt_cox <- proof_genes_pt_cox_recurrence
    test_data <- test_data_recurrence
  } else{
    analysis <- "survival"
    final_fit  <-  final_fit_survival
    train_data <-  train_data_survival
    cox_model  <-  cox_model_survival
    proof_genes_pt_cox <- proof_genes_pt_cox_survival
    test_data <- test_data_survival
  }
  
  
  
  # 1.- Lasso impact analysis -----------------------------------------------
  
  #> First we do a plot that evaluates how genes are penalized on a Lasso model when lambda increases
  #> For that first we extract coefficients from the final fit
  
  
  tidy_coeffs <- tidy(extract_fit_engine(final_fit))
  
  # 1.2 We then do an inverse log so as that to plot the estimates with relationship with the lambda
  
  plot_data <- tidy_coeffs %>%
    mutate(log_lambda = log(lambda))
  
  # 1.3 Filter for just the endpoints (where log_lambda is at its maximum) so that when we assign labels it only assigns to that point
  
  endpoints <- plot_data %>%
    group_by(term) %>%
    filter(log_lambda == min(log_lambda)) %>%
    ungroup()
  
  # 1.4 Plot
  
  lasso_plot <-
    ggplot(plot_data,
           aes(
             x = log_lambda,
             y = estimate,
             group = term,
             color = term
           )) +
    geom_line(alpha = 0.7) + # Make lines slightly transparent so labels pop
    geom_text_repel(
      data = endpoints,
      aes(label = term),
      size = 3,
      hjust = 0,
      nudge_x = -3,
      # Pushes labels further to the right
      direction = "both",
      # Forces labels to only move up/down to avoid overlap
      segment.color = "grey70",
      box.padding = 1,
      # Reduces space around the text boxes
      max.overlaps = Inf         # Ensures every gene gets a name, even if it's crowded
    ) +
    theme(plot.margin = margin(r = 500)) +
    scale_color_viridis_d(option = "turbo") + # High-contrast palette for many genes
    theme_minimal() +
    theme(legend.position = "none") + # Legend is redundant since names are on lines
    expand_limits(x = max(plot_data$log_lambda) - 1.5) # Make room for the names
  
  
  # 2.- Brier score, Schonfeild and martingale residuals --------------------
  
  # 1 Add a row ID so we can find these patients later
  
  id <- train_data %>%
    rownames_to_column("PATIENT_ID") %>%
    mutate(EVENT_STAT = as.numeric(as.character(EVENT_STAT)))
  
  # 1.2 Utilize the fitted object to make predictions based on time
  
  model_diagnostics <- augment(
    final_fit,
    new_data = id,
    eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
  )
  
  
  # 2.1 Brier score
  
  model_diagnostics %>%
    brier_survival(truth = surv_obj, .pred)
  
  
  # 2.2 Martingale and Schofeild residuals
  
  print(cox.zph(cox_model))
  
  print(ggcoxzph(cox.zph(cox_model)))
  
  print(ggcoxdiagnostics(
    cox_model,
    type = "martingale",
    linear.predictions = FALSE,
    ggtheme = theme_bw()
  ))
  
  # 3.- Gene global contribution to score -----------------------------------
  
  # 3.1 Create a data frame that cntains the estimates, if its protectigve or risl and the importance in positive values
  
  coef_df <- tidy(final_fit) %>%
    filter(estimate != 0) %>%
    mutate(direction = ifelse(estimate < 0, "Protective", "Risk"),
           importance = estimate^2)
  
  # 3.2 Plot descending importance
  
  explainability_plot <-
    ggplot(coef_df, aes(
      x = estimate,
      y = reorder(term, importance),
      fill = direction
    )) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c(
      "Risk" = alpha("#b545c43f", 0.8),
      "Protective" = alpha("#8324d2cd", 0.8)
    )) +
    labs(x = "Coefficient (log hazard)", y = NULL, title = stringr::str_to_title(analysis)) +
    theme_minimal() + 
    theme(
      legend.position = ifelse(analysis == "recurrence", "right", "none")
    )
  
  
  
  # 6.- Observe distributions and shapiro -----------------------------------
  
  
  # 6.1 Shapiro wilk test of desired distribution
  
  shapiro.test(train_data$risk_score)
  
  # 6.2 Plot distribution
  
  ggplot(train_data, aes(x = risk_score)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "steelblue",
      color = "black"
    ) +
    geom_density(color = "#81dcff", size = 1) +
    theme_minimal()
  
  
  # 6.3 Plot distribution based on event stat
  
  ggplot(test_data, aes(x = risk_score, fill = factor(
    EVENT_STAT, labels = c("Alive", "Diceased")
  ))) +
    geom_density(alpha = 0.4) +
    theme_linedraw() +
    labs(
      fill = "Event",
      x = "Risk score",
      y = "Density",
      title = "Distribution of score with respect to event"
    )
  
  if (i == 1) {
    lasso_plot_rec <- lasso_plot
    explainability_plot_rec <- explainability_plot

  } else{
    lasso_plot_surv <- lasso_plot
    explainability_plot_surv <- explainability_plot
  }
  

}

(explainability_plot_surv |
explainability_plot_rec)

