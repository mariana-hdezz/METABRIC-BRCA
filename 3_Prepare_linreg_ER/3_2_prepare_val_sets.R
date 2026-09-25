library(dplyr)
library(tidyr)
library(tibble)

proof_genes_surv <- readRDS("./output_data/proof_genes_surv.RDS")
proof_genes_rec  <- readRDS("./output_data/proof_genes_rec.RDS")

# finall object GSE2034 ---------------------------------------------------

gene_expres_matrix_gse2034 <- readRDS("./output_data/gene_expres_matrix_gse2034.RDS")
metadata_gse_2034_er_pos   <- readRDS("./output_data/metadata_gse_2034_er_pos.RDS")

# 1.2 Object with all the patients and expression of only the genes of interest

gene_expres_matrix_gse2034_filt <- gene_expres_matrix_gse2034[,colnames(gene_expres_matrix_gse2034) %in% proof_genes_rec]

# Keep only genes to analyze and join with the metadata_gse_2034_er_pos that was cleaned earlier

proof_genes_pt_gse2034 <- 
  gene_expres_matrix_gse2034_filt[, proof_genes_rec] %>% # Keep only the genes to test
  as.data.frame() %>% 
  rownames_to_column("file_name") %>%
  left_join(metadata_gse_2034_er_pos, by = "file_name") %>% # Join with metadata
  mutate(
    surv_obj = Surv(time = EVENT_MON, event =  EVENT_STAT) # Create the Survival Object inside the dataframe
  ) %>% 
  column_to_rownames("file_name") %>% 
  dplyr::select(- index, # Select only the columns of interest and discard the rest
                - id,
                - ER_STAT,
                - BRAIN_REL,
                -LYMPH,
                - sampleNames
                ) %>% 
  filter(EVENT_MON >= 2 & EVENT_MON <= 180)



proof_genes_pt_gse2034 <- proof_genes_pt_gse2034 %>% 
  mutate(across(-c(EVENT_MON, EVENT_STAT, surv_obj), ~ as.vector(scale(.x))))

saveRDS(proof_genes_pt_gse2034, "./output_data/proof_genes_pt_gse2034.RDS")

# Finall object GSE96058 --------------------------------------------------

counts_data_gse96058_erpos <- readRDS("./output_data/counts_data_gse96058_erpos.RDS")
metadata_gse96058_er_pos <- readRDS("./output_data/metadata_gse96058_er_pos.RDS")

# 4.2 Object with all the patients and expression of only the genes of interest

counts_data_gse96058_erpos <- counts_data_gse96058_erpos[, colnames(counts_data_gse96058_erpos) %in% proof_genes_surv]


# 4.4 Keep only genes that are to be evaluated

counts_data_gse96058_erpos <- counts_data_gse96058_erpos[ , proof_genes_surv]

# 4.5 Object with surv object, all of the genes and the event columns

proof_genes_pt_gse96058 <- 
  counts_data_gse96058_erpos %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata_gse96058_er_pos, by = "title") %>% 
  mutate(surv_obj = Surv(time = EVENT_MON, event = EVENT_STAT, type = "right")) %>% 
  dplyr::select(all_of(proof_genes_surv),
                EVENT_STAT,
                EVENT_MON,
                title,
                surv_obj) %>% 
  column_to_rownames("title")


proof_genes_pt_gse96058 <- 
  proof_genes_pt_gse96058 %>% 
  mutate(across(-c(EVENT_MON, EVENT_STAT, surv_obj), ~ as.vector(scale(.x))))

saveRDS(proof_genes_pt_gse96058, "./output_data/proof_genes_pt_gse96058.RDS")

# Final objectS TCGA ------------------------------------------------------

counts_data_tcga <- readRDS("./output_data/counts_data_tcga.RDS")
refined_data_unique <- readRDS("./output_data/refined_data_unique.RDS")


for (i in c(1:2)) {
  if (i == 1) {
    gene_signature <- proof_genes_surv
    analysis <- "SURVIVAL"
    
  } else{
    gene_signature <- proof_genes_rec
    analysis <- "RECURRENCE"
    
  }
  
  
  # 4.2 Object with all the patients and expression of only the genes of interest
  
  proof_genes_pt_tcga <- counts_data_tcga[rownames(counts_data_tcga) %in% gene_signature, ]
  
  # 5.- Prepare object for validation ---------------------------------------
  
  # 5.1 Transpose first so samples are rows
  
  tcga_transposed <-
    t(proof_genes_pt_tcga)
  
  
  # 5.3 Assign log to the object for ML
  
  proof_genes_pt_tcga <-
    tcga_transposed %>% as.data.frame()
  
  
  # 5.4 Check that the patients are in the same order
  
  refined_data_unique <-
    refined_data_unique[refined_data_unique$sampleID %in% rownames(proof_genes_pt_tcga), ]
  
  all(rownames(proof_genes_pt_tcga) == refined_data_unique$sampleID)
  
  # 5.5 Add a column of EVENT as a binary term for it to be the outcome
  if (analysis == "SURVIVAL") {
    proof_genes_pt_tcga <-
      proof_genes_pt_tcga %>%
      rownames_to_column("sampleID") %>%
      left_join(refined_data_unique, by = "sampleID") %>%  # Join counts with metadata
      column_to_rownames("sampleID") %>%
      dplyr::select(
        all_of(gene_signature),
        # Keep all the genes to ve evaluated and the oucome variables
        SURVIVAL_MON,
        # SURVIVAL_MON for survival and RECURRENCE_MON for recurrence
        SURVIVAL
      ) %>% # SURVIVAL for survival and RECURRENCE for recurrence
      dplyr::rename(
        EVENT_STAT = SURVIVAL,
        # Rename to common term (EVENT_STAT for event and EVENT_MON for time of follow up)
        EVENT_MON = SURVIVAL_MON
      ) %>%
      mutate(EVENT_STAT = as.numeric(EVENT_STAT),
             EVENT_MON = as.numeric(EVENT_MON)) %>%
      as.data.frame() %>%
      mutate(surv_obj =  Surv(
        # Create survival object
        time  = EVENT_MON,
        event = EVENT_STAT
      )) %>%
      filter(EVENT_MON > 0)
    
  } else{
    proof_genes_pt_tcga <-
      proof_genes_pt_tcga %>%
      rownames_to_column("sampleID") %>%
      left_join(refined_data_unique, by = "sampleID") %>%  # Join counts with metadata
      column_to_rownames("sampleID") %>%
      dplyr::select(
        all_of(gene_signature),
        # Keep all the genes to ve evaluated and the oucome variables
        RECURRENCE_MON,
        # SURVIVAL_MON for survival and RECURRENCE_MON for recurrence
        RECURRENCE
      ) %>% # SURVIVAL for survival and RECURRENCE for recurrence
      dplyr::rename(
        EVENT_STAT = RECURRENCE,
        # Rename to common term (EVENT_STAT for event and EVENT_MON for time of follow up)
        EVENT_MON = RECURRENCE_MON
      ) %>%
      mutate(EVENT_STAT = as.numeric(EVENT_STAT),
             EVENT_MON = as.numeric(EVENT_MON)) %>%
      as.data.frame() %>%
      mutate(surv_obj =  Surv(
        # Create survival object
        time  = EVENT_MON,
        event = EVENT_STAT
      )) %>%
      filter(EVENT_MON > 0)
    
  }
  
  
  
  if (analysis == "RECURRENCE") {
    proof_genes_pt_tcga <-
      proof_genes_pt_tcga %>%
      filter(EVENT_MON >= 2 & EVENT_MON <= 180)
    print("Recurrene singature, filtered")
  } else{
    print("Survival signature, no filter")
  }
  
  
  proof_genes_pt_tcga <- proof_genes_pt_tcga %>%
    mutate(across(-c(EVENT_MON, EVENT_STAT, surv_obj), ~ as.vector(scale(.x))))
  
  saveRDS(
    proof_genes_pt_tcga,
    paste0(
      "./output_data/",
      "proof_genes_pt_tcga_",
      tolower(analysis),
      ".RDS"
    )
  )

  
}
