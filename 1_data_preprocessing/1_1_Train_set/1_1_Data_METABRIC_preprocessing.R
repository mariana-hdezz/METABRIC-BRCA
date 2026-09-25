library(dplyr)
library(tidyr)
library(tibble)


# 0.- Data downloading ----------------------------------------------------

directory <- readline("Directory for download of METABRIC data")

options(timeout = 600)

download.file(
  url = "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz", 
  destfile = paste0(directory, "brca_metabric.tar.gz")
  )

untar(paste0(directory, "brca_metabric.tar.gz"), exdir = directory)

# 1.- Data loading --------------------------------------------------

microarray_data <- read.table(
  paste0(
    paste0(directory, "brca_metabric/"), 
  "data_mrna_illumina_microarray.txt"
  ),
                              header = TRUE)

# 1.1.2 Metadata

clinical_data <- read.delim(
  paste0(
    paste0(directory, "brca_metabric/"), 
  "data_clinical_patient.txt"
  ),
  comment.char = "#",
  stringsAsFactors = FALSE)


# 1.1.3 Check if rownames in clinical data have the same nomenclature as the colname sin microarray_data, if they differ by a point or dash (MB.0000 vs MB-0000) use

clinical_data$PATIENT_ID <- gsub("-", ".", clinical_data$PATIENT_ID)

# 1.1.4 Object without ENTREZ and HUGO columns

microarray_data_num <- microarray_data[,-(1:2)] 

# 1.2.- Eliminate duplicates and Hugo IDs -----------------------

# 1.2.1 There are Hugo IDs duplicated, and thus we search for variance and keep the duplicate with higher variance

# 1.2.2 Variance

microarray_data.var <- apply(microarray_data_num, 1, var)

# 1.2.3 Asign new column with variance

microarray_data$variance <- microarray_data.var

# 1.2.4 Only maintain the version of the gene duplicate with higher variance

microarray_data_unique <- microarray_data %>% # Initial data
  group_by(Hugo_Symbol) %>% # Group by HUGO symbol
  slice_max(order_by = variance, n = 1) %>% # Order by variance and keep the highest
  ungroup()


# 1.2.5 Convert from tibble to data frame

microarray_data_unique <- as.data.frame(microarray_data_unique)

# 1.2.6 Asign Hugo symbol as rownames

rownames(microarray_data_unique) <- microarray_data_unique$Hugo_Symbol

# 1.2.7 Delete the columns with hugo, entrez and variance.

microarray_data_unique$Hugo_Symbol <- NULL
microarray_data_unique$Entrez_Gene_Id <- NULL
microarray_data$variance <- NULL

counts_data <- microarray_data_unique
counts_data$variance <- NULL

##############################################

# 2.- METADATA -----------------------


# 2.1 Metadata that contains only the patients of whom we have RNA counts

metadata <- clinical_data[clinical_data$PATIENT_ID %in% colnames(microarray_data_num),]


metadata  <-   metadata %>% 
  mutate(ER_POS = ifelse(ER_IHC == "Positve", # Create ER+ column with 0 and 1 depending on status
                         1,
                         0),  
         SURVIVAL_STAT = as.numeric(gsub(":.*", "", metadata$OS_STATUS)), # 2.2 Add column of survival state where 1 is dead and 0 is alive
         RECURR_STAT = as.numeric(gsub(":.*", "", metadata$RFS_STATUS))) %>% 
  mutate(HER2_POS = ifelse(THREEGENE == "HER2+",
                           1,
                           0),
         HER2_POS = as.numeric(HER2_POS),
         HER2 = HER2_SNP6, 
         LYMPH = LYMPH_NODES_EXAMINED_POSITIVE,
         AGE = as.numeric(AGE_AT_DIAGNOSIS ),
         MENO = INFERRED_MENOPAUSAL_STATE    ,
         HORMONE = HORMONE_THERAPY ,
         CHEMO = CHEMOTHERAPY,
         SURGERY = BREAST_SURGERY,
         PAM50 = CLAUDIN_SUBTYPE,
         NPI = NPI,
         HIST = HISTOLOGICAL_SUBTYPE,
         RADIO = RADIO_THERAPY,
  )



# 2.4 Object of patients that are alive or  diseased by breast cancer

alive_brca_death <- metadata %>% 
  filter(VITAL_STATUS != "Died of Other Causes") # Relevant for recurrence



# 2.5 Complete metadata only for ER+ patients

metadata_ER_POS <- metadata %>% 
  as.data.frame() %>% 
  filter(ER_IHC == "Positve")  # Keep only ER+ patients

# 2.5.1 Complete metadata only for ER+ useful for recurrence (min 2 months, max 180 months)

metadata_ER_POS_REC <- metadata_ER_POS %>% 
  as.data.frame() %>% 
  mutate(EVENT_STAT = as.numeric(RECURR_STAT),
         EVENT_MON = as.numeric(RFS_MONTHS)) %>% 
  filter(EVENT_MON <= 180 & EVENT_MON >= 2) 

# 2.5.2 Complete metadata only for ER+ patients who DIED OF BRCA or are alive (useful for survival)

metadata_ER_POS_SURV <- metadata_ER_POS %>% 
  as.data.frame() %>% 
  filter(VITAL_STATUS != "Died of Other Causes") %>% # Only remove the ones who died of other causes
  mutate(EVENT_STAT = as.numeric(SURVIVAL_STAT),
         EVENT_MON = as.numeric(OS_MONTHS)) 




if(dir.exists("./output_data/")){
  "Output adata directory already exists"
}else{
  dir.create("./output_data/")
}


if(dir.exists("./results/")){
  "results directory already exists"
}else{
  dir.create("./results/")
} 


saveRDS(counts_data, "./output_data/counts_data.RDS")
saveRDS(metadata_ER_POS_REC, "./output_data/metadata_ER_POS_REC.RDS")
saveRDS(metadata_ER_POS_SURV, "./output_data/metadata_ER_POS_SURV.RDS")


rm(list = ls())
gc()