library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(survival)
library(paletteer)
library(UCSCXenaTools)
library(ggplot2)


# 1.- Loading data --------------------------------------------------------

directory <- readline("Directory for download of TCGA data")

# 1.1 Query con base a RNA-Seq y a 3 casos y 3 controles

tcga_rna <- GDCquery(
  "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  access = "open",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"
)

GDCdownload(
  tcga_rna,
  method = "api",
  files.per.chunk = 5,
  directory = directory
)

# 1.2 Prepare data for usage

tcga_brca_data <- GDCprepare(tcga_rna, directory = directory)

# 1.3 Count matrix

brca_matrix <- assay(tcga_brca_data, "fpkm_unstrand")

genes <- rowData(tcga_brca_data)[colnames(rowData(tcga_brca_data)) %in% c("gene_name")]

rm("tcga_brca_data")

# Assign genes as names

brca_matrix <- merge(brca_matrix, genes, by = 0)

brca_matrix <- brca_matrix %>% 
    dplyr::select(- Row.names)

# 5.2 Log2 Transform

brca_matrix_log <- as.data.frame(brca_matrix)

numeric_cols <- sapply(brca_matrix_log, is.numeric)

brca_matrix_log[, numeric_cols] <- sapply(brca_matrix_log[, numeric_cols], function(x) log2(x + 1))

brca_data <- brca_matrix_log

# 1.4.2 Convert the dots to slashes

colnames(brca_data) <- gsub("\\.", "-", colnames(brca_data))


# 1.5.- Eliminate duplicates ----------------------------------------------


# 1.5 Extract sample type from TCGA barcode

# 1.5.2 Select the 14th - 16th value which correspond to sample type codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# this to then select only the samples that correspond to primary tumors

sample_type_full <- substr(colnames(brca_data), 14, 16)

# 1.5.3 Maintain only primary tumor samples so as to avoid duplicates

# 1.5.3.1 Keep only counts that correspond to primary tumor

brca_data <- brca_data[, sample_type_full == "01A" | colnames(brca_data) %in% "gene_name"]

# 1.5.4 Assign to new object that will be modified to eliminate duplicates

brca_data2 <- brca_data

# 1.5.4.2 Keep the names up until the -01 so as to have it in the same nomenclature as the metadata

colnames(brca_data2) <- substr(colnames(brca_data2), 1, 15)

# 1.5.4.3 We can see that there are 5 patients with 2 samples of the same tumor

names(brca_data)[substr(colnames(brca_data), 1, 15) %in% names(brca_data2)[duplicated(names(brca_data2))]]

# 1.5.5 So we keep only the patient sample with highest variance

# 1.5.5.2 Calculate variance

sample_variance <- apply(brca_data2[,!colnames(brca_data2) %in% "gene_name"], 2, var, na.rm = TRUE)

# 1.5.5.3 Match variance with the names

df <- data.frame(
  sample = colnames(brca_data2[,!colnames(brca_data2) %in% "gene_name"]),
  equivalent = colnames(brca_data[,!colnames(brca_data) %in% "gene_name"]),
  variance = sample_variance # column with the variance
)

# 1.5.5.3 Group by and keep only the sample with highest variance

selected_samples <- df %>%
  group_by(sample) %>% # Group by the reduced name (TCGA-XX-XXXX-01)
  slice_max(order_by = variance,
            n = 1,
            with_ties = FALSE) %>% # Keep only the highest variance
  pull(equivalent) # Extract the complete name (TCGA-XX-XXXX-01A-XXX-XXXX-XX)

# 1.6 Keep only the counts of the unique

brca_data_unique <- brca_data[, c(selected_samples, "gene_name")]

# 1.7 Convert back to names compatible with metadata

colnames(brca_data_unique) <- substr(colnames(brca_data_unique), 1, 15)

# 2.- Metadata ------------------------------------------------------------

# 2.1 Generate and Query

data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>%
  XenaQuery()

# 2.2 Download

xe_download <- XenaDownload(data_query, destdir = "D:/tcga/GDCdata/Metadata")

# 2.3 Prepare (Load) the data

brca_clinical <- XenaPrepare(xe_download)

# 2.4 Create the Recurrence variables


refined_data <-
  brca_clinical %>%
  mutate(
    SURVIVAL = ifelse(vital_status == "DECEASED", 1, 0),
    SURVIVAL_MON = ifelse(
      !is.na(days_to_death) & days_to_death >= 0,
      days_to_death,
      # If they died they have this parameter
      days_to_last_followup # Else they have this parameter
    ) / 30.4166667,
    RECURRENCE = ifelse(
      new_neoplasm_event_type %in% c("Locoregional Recurrence", "Distant Metastasis"),
      1,
      0
    ),
    RECURRENCE_MON = ifelse(
      RECURRENCE == 1,
      days_to_new_tumor_event_after_initial_treatment,
      coalesce(days_to_last_followup, days_to_death)
    ) / 30.4166667
  ) %>%
  filter(
    breast_carcinoma_estrogen_receptor_status == "Positive" # Only ER+,
  ) %>%
  mutate(
    HER2 = lab_proc_her2_neu_immunohistochemistry_receptor_status,
    LYMPH = as.numeric(lymph_node_examined_count),
    PAM50 = PAM50Call_RNAseq,
    AGE = as.numeric(Age_at_Initial_Pathologic_Diagnosis_nature2012),
    RADIO = radiation_therapy,
    SURGERY = factor(breast_carcinoma_primary_surgical_procedure_name),
    NEO = history_of_neoadjuvant_treatment,
    OTHER_TX = additional_pharmaceutical_therapy,
    TARG_TX = targeted_molecular_therapy,
    HER2 = HER2_Final_Status_nature2012,
    MENO = menopause_status,
    HIST = histological_type,
    INTCLUST = Integrated_Clusters_with_PAM50__nature2012
  )


# 2.5 To maintain the samples of the primary tumor and eliminate duplicates as done with counts data

sample_type_full2 <- substr(refined_data$sampleID, 14, 15)

refined_data_unique <- refined_data

refined_data_unique <- refined_data_unique[sample_type_full2 == "01", ]


# 2.6 Count object that corresponds to the metadata patients

brca_data_filtered <- brca_data_unique[, colnames(brca_data_unique) %in% c(refined_data_unique$sampleID, "gene_name")]

# 2.7 Since there are missing patients since the first data we also filter out patients in metadata that arent on counts

refined_data_unique <- refined_data_unique[refined_data_unique$sampleID %in% colnames(brca_data_filtered), ]



# 3.- Deleting duplicates and asigning ensembl as rownames ----------------

# 3.2 Variance

numeric_data <- brca_data_filtered %>%
  dplyr::select(where(is.numeric))

brca_data_filtered$variance <- apply(numeric_data, 1, var)


# 3.2.2 Only mantain the version of the gene duplicate with higher variance

counts_data_tcga <- brca_data_filtered %>% # Initial data
  group_by(gene_name) %>% # Group by ensembl
  slice_max(order_by = variance,
            n = 1,
            with_ties = FALSE) %>% # Order by variance and keep the highest
  ungroup() %>%
  filter(!variance == 0) %>%
  dplyr::select(-variance) %>% # Delete variance, and both of the ensembl ids columns
  filter(!is.na(gene_name)) %>%  # Delete those that had a NA in symbol
  column_to_rownames("gene_name")


saveRDS(counts_data_tcga, "./output_data/counts_data_tcga.RDS")
saveRDS(refined_data_unique,
        "./output_data/refined_data_unique.RDS")


rm(list = ls())
gc()
