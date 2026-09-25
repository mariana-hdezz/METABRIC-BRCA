library(limma)
library(tidyverse)


outlier_cause_surv <- readRDS("./output_data/outlier_cause_surv.RDS")
outlier_cause_rec <- readRDS("./output_data/outlier_cause_rec.RDS")

counts_data <- readRDS("./output_data/counts_data.RDS")

for (i in c(1, 2)) {
  if (t == 1) {
    metadata  <- outlier_cause_rec
    analysis <- "recurrence"
    
  } else{
    metadata  <-  outlier_cause_surv
    analysis <- "survival"
  }
  
# 1.- Dividing by BIAS nodes ------------------------------

# 1.1 Generate column corresponding to BIAS nodes, those that have 0 in one group and those with more than 0 in another

  for (t in c(1, 2)) {
  
col_data <- metadata %>%
  mutate(BIAS = quadrant) %>% 
  filter(BIAS == t | BIAS == 3) %>%
  mutate(BIAS = ifelse(BIAS == 3, 0, 1)) %>% 
  dplyr::select(PATIENT_ID, BIAS) %>% # Create only the object to use for Limma
  column_to_rownames("PATIENT_ID")

if(t == 1){
  quadrant <- "unexcpected"
}else{
  quadrant <- "exceptional"
}

# 2.- Differential expression -----------------------------------------------


# 2.2 Data counts of the patients that had BIAS node information in the metadata

count_data <- counts_data[colnames(counts_data) %in% rownames(col_data)]

# 2.2.2 Making shure they are in the same order

count_data <- count_data[match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))


# 2.3 Generate limma object

library(limma)

# 2.4 Design based on object separating on BIAS nodes

col_data$BIAS <- as.factor(col_data$BIAS)

design <- model.matrix(~ 0 + BIAS, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design)) 

# 2.5 Fit

fit <- lmFit(count_data, design)

# 2.5.2 Contrast matrix comparing BIAS 0 to > 0 BIAS

contrast.matrix <- makeContrasts(BIAS1 - BIAS0,
                                 levels = design)

# 2.5.3 Fit based on contrasts

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit)

# 2.6 Results

res <- topTable(fit, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) # 0.1
res_sig

intersect(proof_genes, rownames(res_sig))

source("./6_bias_and_misclassification/6_3_3_gsea_missclass.R")

write.csv(res, paste0("./results/", "diff_expr_missclass_", quadrant, "_", analysis, ".RDS"))
write.csv(gse_df, paste0("./results/", "gsea_missclass_", quadrant, "_", analysis, ".RDS"))

  }
}
