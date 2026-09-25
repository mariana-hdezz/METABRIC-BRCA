

# 2.- Identifying cause of misclassification ------------------------------

outlier_cause_surv <- readRDS("./output_data/outlier_cause_surv.RDS")
outlier_cause_rec <- readRDS("./output_data/outlier_cause_rec.RDS")

for (i in c(1, 2)) {
  if (i == 1) {
    analysis <- "recurrence"
    outlier_cause  <- outlier_cause_rec
    
  } else{
    analysis <- "survival"
    outlier_cause  <-  outlier_cause_surv
  }
  
  print(analysis)
  
  # 2.1 Comparison between groups and p val adjustment
  
  for (a in c(1, 2)) {
    # 2.2 Subset data depending on cuadrant analysed
    
    df_misclassification <- subset(outlier_cause, quadrant %in% c(a, 3)) # First Q1 vs Q3 (unexpected vs correct) thenb Q2 vs Q3 (exceptional vs correct)
    
    p_multiple_eval_misclass <- c(
      intclust = chisq.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$INTCLUST
        ),
        simulate.p.value = TRUE
      )$p.value,
      # Many variables unfit for fishers but because of possible sparsity we p simulate
      histology = chisq.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$HISTOLOGICAL_SUBTYPE
        ),
        simulate.p.value = TRUE
      )$p.value,
      # Same
      hormone = fisher.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$HORMONE_THERAPY
        )
      )$p.value,
      intcluster  = fisher.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$intcluster
        )
      )$p.value,
      # Binary withprobaiblity of <5 variables in cells
      chemo    = fisher.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$CHEMOTHERAPY
        )
      )$p.value,
      radio    = fisher.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$RADIO_THERAPY
        )
      )$p.value,
      claudin  = chisq.test(
        table(
          df_misclassification$quadrant,
          df_misclassification$CLAUDIN_SUBTYPE
        ),
        simulate.p.value = TRUE
      )$p.value,
      nodes    = kruskal.test(LYMPH_NODES_EXAMINED_POSITIVE ~ quadrant, data = df_misclassification)$p.value,
      # Continuous non parametric data
      age      = kruskal.test(AGE_AT_DIAGNOSIS ~ quadrant, data = df_misclassification)$p.value,
      NPI      = kruskal.test(NPI ~ quadrant, data = df_misclassification)$p.value
    )
    cat("\n\n\n")
    # 2.3 Adjust p value with bonferroni holm and print which groups are being comopared
    
    adj_p_misclass_event <- p.adjust(p_multiple_eval_misclass, method = "holm")
    print(paste(analysis, " Q", a, " vs Q3 (Standard) Holm Adjusted P-values"))
    options(scipen = 999)
    print(round(adj_p_misclass_event, 5))
  

    
    
    for (e in c("CLAUDIN_SUBTYPE", "INTCLUST")) {
      # 2.4 Table with characteristic to evlauate its residuals
      
      misclass_table <- table(df_misclassification$quadrant, df_misclassification[[e]])
      misclass_table <- misclass_table[drop = TRUE]
      
      # 2.5 Run chisq with simulated p values
      
      chisq_obj_exp <- chisq.test(misclass_table, simulate.p.value = TRUE)
      
      # 2.6 Obtain residuals
      print(chisq_obj_exp)
      print(chisq_obj_exp$residuals)
    }
    
  }
  
  
  # 3.- Regressions ---------------------------------------------------------
  
  # 3.1 Firth penalized logistic regression on unexpectedevent group
  
  outlier_cause$PAM50 <- relevel(as.factor(outlier_cause$PAM50), "LumB")
  
  clean_unexpected_model_firth <- logistf(
    unexpected_event ~ HORMONE_THERAPY + intcluster + PAM50 + CHEMOTHERAPY + RADIO_THERAPY + AGE_AT_DIAGNOSIS + NPI ,
    data = outlier_cause
  )
  print(summary(clean_unexpected_model_firth))
  # 3.1.2 OR
  print(exp(coef(clean_unexpected_model_firth)))
  # 3.1.3 CI OR
  print(exp(confint(clean_unexpected_model_firth)))
  
  # 3.2 MFirth penalized regression on exceptional group
  
  
  clean_exceptional_model_firth <- logistf(
    exceptional ~ HORMONE_THERAPY + intcluster + PAM50 + CHEMOTHERAPY + RADIO_THERAPY + AGE_AT_DIAGNOSIS + NPI ,
    data = outlier_cause
  )
  print(summary(clean_exceptional_model_firth))
  # 3.2.2 OR
  print(exp(coef(clean_exceptional_model_firth)))
  # 3.2.3 OR with CI
  print(exp(confint(clean_exceptional_model_firth)))
  
  cat("\n\n\n\n\n\n\n\n\n")
}
