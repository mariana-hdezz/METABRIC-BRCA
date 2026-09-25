library(tidymodels)
library(censored)
library(survminer)
library(broom)
library(survival)
library(timeROC)
library(flextable)

# Load

train_data_survival <- readRDS("./output_data/train_data_survival.RDS")
train_data_recurrence <- readRDS("./output_data/train_data_recurrence.RDS")
best_params_surv <- readRDS("./output_data/best_params_survival.RDS")
best_params_recu <- readRDS("./output_data/best_params_recurrence.RDS")

for (i in c(1:2)) {
  if (i == 1) {
    train_data <- train_data_survival
    analysis <- "survival"

    
  } else{
    train_data <- train_data_recurrence
    analysis <- "recurrence"

  }
  
  
  # 2.- Preparing recipe and model ------------------------------------------
  
  # 2.1 Recipe
  
  lr_rec <- recipe(surv_obj ~ ., data = train_data) %>%
    update_role(EVENT_MON, EVENT_STAT, new_role = "non_predictor") %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>% # Eliminates variables with a single value
    step_nzv(all_predictors()) # Eliminates highly sparse variables
  
  
  # 2.2 Model
  
  lr_mod <- proportional_hazards(
    penalty = tune(),
    # lambda establishes the severith of the penalty
    mixture = tune()     # alpha establishes the type, 1 being lasso, 0 being ridge, and 0.5 being elasticnet
  ) %>%
    set_engine("glmnet", cox.ties = "breslow") # Engine that permits penalizing by elasticnet, ridge, and lasso
  
  # 2.3 Workflow
  
  lr_wf <- workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(lr_rec)
  

    set.seed(123)
    
    # 3.1 Folds for evaluating with resamples
    
    folds <- vfold_cv(
      train_data,
      v = 10,
      strata = EVENT_STAT
    )
    
    # 3.2 Grid for penalizing range
    
    grid <- grid_regular(
      penalty(range = c( - 4, 1)),   
      mixture(range = c(0, 1)),
      levels = 10
    )
    
    # 3.3 Running the different penalization methods
    
    res_ml <- tune_grid(
      lr_wf,
      resamples = folds,
      grid = grid,
      metrics = metric_set(concordance_survival), # Evaluates the different penalizing methods by c score
      control = control_grid(save_pred = TRUE)
    )
    
    # 3.3.1 Observe metrics
    
    collect_metrics(res_ml)
    
    # 3.3.2 Object with best parameters for penalizing
    
    best_params <- select_best(res_ml, metric = "concordance_survival")
    print(best_params)
  
  # 4.- Actual training -----------------------------------------------------
  
  # 4.1 Final workflow with the selecting the best parameter tested previously
  
  final_wf <- finalize_workflow(lr_wf, best_params)
  
  # 4.2 Final fit with training data
  
  final_fit <- fit(final_wf, data = train_data)
  
  # 4.2.2 Observing genes that are maintained after penalziation
  
  coef_tbl <- tidy(final_fit) %>%
    filter(estimate != 0) %>%
    arrange(desc(abs(estimate)))
  
  cat(coef_tbl$term, sep = ", ")
  
  # 4.3 Predictions on train data
  
  train_pred <- predict(final_fit, new_data = train_data, type = "linear_pred")
  
  train_pred <-
    train_pred %>%
    as.data.frame()
  
  # 4.3.2 Object with train data and its predicted scores
  
  train_data2 <-
    train_data %>%
    mutate(risk_score = train_pred$.pred_linear_pred)
  
  # 4.4 Using those scores to calculate the cutpoint
  
  true_cut <-
    surv_cutpoint(
      data = train_data2,
      time = "EVENT_MON",
      event = "EVENT_STAT",
      variables = "risk_score"
    )
  
  
  saveRDS(train_data2, paste0("./output_data/train_data_", analysis, ".RDS"))
  saveRDS(final_fit, paste0("./output_data/final_fit_", analysis, ".RDS"))
  saveRDS(true_cut, paste0("./output_data/true_cut_", analysis, ".RDS"))
  
}
