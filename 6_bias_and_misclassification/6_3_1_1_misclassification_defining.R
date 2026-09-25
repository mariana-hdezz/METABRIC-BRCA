library(dplyr)
library(tidyr)
library(tibble)
library(tidymodels)
library(survival)
library(logistf)



final_fit_survival    <- readRDS("./output_data/final_fit_survival.RDS")
train_data_recurrence <- readRDS("./output_data/train_data_recurrence.RDS")

final_fit_recurrence  <- readRDS("./output_data/final_fit_recurrence.RDS")
train_data_survival   <- readRDS("./output_data/train_data_survival.RDS")

metadata_ER_POS_REC <- readRDS("./output_data/metadata_ER_POS_REC.RDS")
metadata_ER_POS_SURV <- readRDS("./output_data/metadata_ER_POS_SURV.RDS")


for (t in c(1:2)) {
  if (t == 1) {
    analysis <- "recurrence"
    final_fit  <- final_fit_recurrence
    train_data <- train_data_recurrence
    ml_metadata <- metadata_ER_POS_REC
  } else{
    analysis <- "survival"
    final_fit  <-  final_fit_survival
    train_data <-  train_data_survival
    ml_metadata <- metadata_ER_POS_SURV
  }
  
  
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
  
  # 1.2.2 Create list to then add the results
  
  outlier_time_list <- list()
  
  no_img <- 0
  
  #> 1.3 For loop that at each desired time point calculates the misclassification scores, extracts patients with high misclassification scores
  #> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time
  
  for (i in c(36, 60, 120)) {
    eval_time <- i
    
    # 1.3 Unnest the predictions and find the biggest outliers on a set point and event
    
    outliers <- model_diagnostics %>%
      dplyr::select(PATIENT_ID, EVENT_MON, EVENT_STAT, .pred, .pred_time) %>%
      unnest(.pred) %>%
      filter(.eval_time == eval_time) %>%
      arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive
    
    
    # 1.4 Object with metadata and score characteristics
    
    outlier_summary <- outliers %>%
      inner_join(ml_metadata,
                 by = "PATIENT_ID",
                 suffix = c("", ".drop")) %>%
      dplyr::select(
        PATIENT_ID,
        .pred_survival,
        EVENT_STAT,
        EVENT_MON,
        CLAUDIN_SUBTYPE,
        LYMPH_NODES_EXAMINED_POSITIVE,
        THREEGENE,
        INTCLUST,
        .eval_time,
        NPI,
        CELLULARITY,
        HISTOLOGICAL_SUBTYPE,
        .pred_time,
        CHEMOTHERAPY,
        RADIO_THERAPY,
        HORMONE_THERAPY,
        BREAST_SURGERY,
        AGE_AT_DIAGNOSIS,
        OS_STATUS,
        OS_MONTHS,
        RFS_MONTHS,
        RFS_STATUS
      )
    
    
    # 1.5 Create misclassification score for defined time
    
    outliers_misclassification <- outlier_summary %>%
      mutate(misclassification_score = (((1 - EVENT_STAT) - .pred_survival)^2) * (EVENT_MON - .eval_time) * (1 - (2 * EVENT_STAT))) %>%
      arrange(desc(misclassification_score))
    
    # 1.6 Identify the highest misclassification patients
    
    q3 <- quantile(outliers_misclassification$misclassification_score,
                   0.75)
    iqr_val <- IQR(outliers_misclassification$misclassification_score)
    upper_fence <- q3 + (1.5 * iqr_val)
    
    # 1.6.1 Obtain mean and sd and then filter baed on patients higher than determined SD
    
    extreme_outliers <-
      outliers_misclassification %>%
      mutate(
        mean_misclassification = mean(misclassification_score),
        sd_misclassification = sd(misclassification_score)
      ) %>%
      filter(misclassification_score > upper_fence)
    
    # 1.6.2 Obtain their IDs
    
    top_misclassification_ids <- extreme_outliers$PATIENT_ID
    
    # 1.7.1 Add a column identifying patients as top misclassification or not
    
    outliers <-
      outliers %>%
      mutate(
        quadrant = case_when(
          PATIENT_ID %in% top_misclassification_ids & EVENT_STAT == 0 ~ 2,
          PATIENT_ID %in% top_misclassification_ids &
            EVENT_STAT == 1 ~ 1,
          TRUE ~ 0
        )
      )
    
    # 1.7.2 Plot
    if (no_img == 0 & eval_time == 60) {
      theme_embedded <- theme_classic(base_size = 15) +
        theme(
          legend.position = c(0.88, 0.2),
          # Adjust coordinates (x, y) from 0 to 1
          legend.background = element_rect(fill = alpha("white", 0.5))
        )
      
      # 1.7.2.1 Plot colored by EVENT_STAT
      
      p1 <- ggplot(outliers,
                   aes(
                     x = EVENT_MON,
                     y = .pred_survival,
                     color = factor(EVENT_STAT),
                     shape = factor(EVENT_STAT)
                   )) +
        geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
        stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
        geom_vline(xintercept = eval_time,
                   linetype = "dashed",
                   color = "red") +
        scale_color_manual(
          values = c("#DB6D00FF", "#920000FF"),
          labels = c(
            ifelse(analysis == "recurrence", "Non-recurrent", "Alive"),
            ifelse(analysis == "recurrence", "Recurrent", "Deceased")
          )
        ) +
        scale_shape_manual(values = c(16, 17),
                           labels = c(
                             ifelse(analysis == "recurrence", "Non-recurrent", "Alive"),
                             ifelse(analysis == "recurrence", "Recurrent", "Deceased")
                           )) +
        labs(
          title = "Event Status Distribution",
          x = paste0("Actual Event Time (Months)", eval_time),
          y = paste0("Predicted ", analysis),
          color = "Event Stat",
          shape = "Event Stat"
        ) +
        theme_embedded
      
      # 1.7.2.2 Plot colored by quadrant
      print(analysis)
      
      p2 <- ggplot(outliers,
                   aes(
                     x = EVENT_MON,
                     y = .pred_survival,
                     color = factor(quadrant),
                     shape = factor(EVENT_STAT)
                   )) +
        geom_point(size = 2, alpha = 0.7) +
        stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) +
        geom_vline(xintercept = eval_time,
                   linetype = "dashed",
                   color = "red") +
        scale_color_manual(
          values = c("#006DDBFF", "#490092FF", "#B66DFFFF"),
          labels = c(
            "Correct",
            "Unexcpected event group",
            "Exceptional outcome group"
          )
        ) +
        scale_shape_manual(values = c(16, 17),
                           labels = c(
                             ifelse(analysis == "recurrence", yes = "Non-recurrent", no = "Alive"),
                             ifelse(analysis == "recurrence", "Recurrent", "Deceased")
                           )) +
        labs(
          title = "Misclassification group Analysis",
          x = paste0("Actual Event Time (Months)", eval_time),
          y = "Predicted Survival",
          color = "Groups",
          shape = "Event Stat"
        ) +
        theme_embedded
      
      
      
      # 1.7.3 Combine and stack
      
    } else{
      
    }
    
    outlier_time_list[[i]] <- top_misclassification_ids
    
  }
  
  
  # 1.8 Obtain patients found on all of the iterations of the for loop as top misclassification patients
  
  misclassification_interesct <- intersect(intersect(outlier_time_list[[36]], outlier_time_list[[60]]) ,
                                           outlier_time_list[[120]])
  
  # 1.8.2 Similar but all unique so even if they appear once we register them
  
  misclassification_diff_id <- unique(c(outlier_time_list[[36]], outlier_time_list[[60]], outlier_time_list[[120]]))
  
  
  # 2.1 Mutate characteristics to evaluate in the further analysis
  
  
  outlier_cause <- id %>%
    inner_join(ml_metadata, by = "PATIENT_ID", suffix = c("", ".drop")) %>%
    mutate(
      # Binarize IntClust to low and high grade
      intcluster = ifelse(INTCLUST %in% c("3", "4ER+", "7", "8"), "Low", "High"),
      quadrant = case_when(
        PATIENT_ID %in% misclassification_diff_id & EVENT_STAT == 1 ~ 1,
        PATIENT_ID %in% misclassification_diff_id &
          EVENT_STAT == 0 ~ 2,
        TRUE ~ 3
      ),
      unexpected_event = ifelse(quadrant == 1, 1, 0),
      exceptional      = ifelse(quadrant == 2, 1, 0)
    )
  
  
  if (analysis == "recurrence") {
    print(analysis)
    saveRDS(outlier_cause, "./output_data/outlier_cause_rec.RDS")
    p1_rec <- p1
    p2_rec <- p2
  } else{
    print(analysis)
    saveRDS(outlier_cause, "./output_data/outlier_cause_surv.RDS")
    p1_surv <- p1
    p2_surv <- p2
  }
  
}


(p1_surv | p2_surv) /
(p1_rec | p2_rec)

