# This script is to perform signature validation on GSE 2034

library(oligo)
library(GEOquery)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(hgu133a.db)
library(phenoTest)
library(clipr)


# 1.- Load data and metadata_gse_2034  ----------------------------------------------

directory <- readline("Directory for download of GSE2034 data")

# 1.1 Get supplementary files

getGEOSuppFiles("GSE2034", baseDir = directory)

# 1.1.2 Untar files
untar(paste0(directory, "GSE2034/GSE2034_RAW.tar"), exdir = paste0(directory, "GSE2034/"))

# 1.1.3 Listing .cel files

cel_files <- list.celfiles(paste0(directory, "GSE2034/"), full.names = TRUE, listGzipped = TRUE)

# 1.1.4 Reading in cel files

pre_raw_data <- read.celfiles(cel_files)

raw_data <- pre_raw_data

# 1.2 metadata_GSE2034 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSE2034&id=40089&db=GeoDb_blob26

pre_metadata <- read_clip_tbl()
metadata_gse_2034 <- pre_metadata

# Object that specifies which GSE is in use

gse_obj <- "GSE2034"


# 3.- Preprocessing metadata_gse_2034  --------------------------------------------------

#3.1 Clean metadata_gse_2034 

metadata_gse_2034  <- metadata_gse_2034  %>%
  mutate(
    EVENT_STAT = ifelse(relapse..1.True. == 1, 1, 0), # Object to evaluate event
    EVENT_STAT = as.numeric( EVENT_STAT),
    EVENT_MON = as.numeric(time.to.relapse.or.last.follow.up..months.), # Object to evaluate time to event
    id = GEO.asscession.number,
    ER_STAT = ER.Status,
    BRAIN_REL = Brain.relapses..1.yes..0.no.,
    LYMPH = lymph.node.status
  ) %>% 
  dplyr::select(EVENT_STAT,
                EVENT_MON,
                ER_STAT,
                id,
                BRAIN_REL,
                LYMPH
                
  )


# 3.2 Object with the names of each file

id <- sampleNames(raw_data)
pData(raw_data)$sampleNames <- id
sample_names <- gsub("\\..*", "", id) # Taking out the .CEL.gz

# 3.2.2 Add the object with the general names to the phenotype data

pData(raw_data)$id <- sample_names

# 3.4 Join metadata_gse_2034  

pData(raw_data) <- 
  pData(raw_data) %>% 
  rownames_to_column("file_name") %>% # To have the full names
  full_join(metadata_gse_2034 , by = "id", keep = FALSE) %>% # Join by names without .CEL.gz
  mutate(comp_file_name = file_name) %>% # Create new column to then add to rownames
  mutate(file_name = comp_file_name) %>% 
  column_to_rownames("comp_file_name") 


metadata_gse_2034  <- pData(raw_data) # Make shure they have the same data this so that metadata also has the identifiers with .CEL.gz to match with the counts



# 4.- Preprocess data -----------------------------------------------------

# Normalize

norm_data <- rma(raw_data) 

# Boxplot after normalization

boxplot(exprs(norm_data), 
        las = 2, 
        main = paste0("RMA normalized - ", gse_obj))

# 4.3 Expression matrix

expr_matrix <- exprs(norm_data)

# 5.- Probe id to symbol --------------------------------------------------

# 5.1 Get mapping to change probe ids to gene names

probe_gene <- AnnotationDbi::select(
  hgu133a.db, # The probe identifiers for this affymetrix
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

# 5.2 Delete NA symbols and duplicate symbols and add to rownames

# 5.2.1 Transform matrix to a tidy data frame and add Symbols

gene_expres_matrix <- 
  expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  inner_join(probe_gene, by = "PROBEID") %>%  # Join with the mapping object
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%   # Remove NAs and empty symbols
  mutate(variance = apply(dplyr::select(., -PROBEID, -SYMBOL), 1, var)) %>%   # 5.2.2 Calculate variance for each probe across all samples
  group_by(SYMBOL) %>% 
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # 5.2.3 Keep only the probe with the highest variance per Gene Symbol
  ungroup() %>%
  dplyr::select(-PROBEID, -variance) %>%  # 5.2.4 Remove columns, add symbol to rownames and reformat to matrix
  column_to_rownames("SYMBOL") %>%
  as.matrix()


# 5.3.3 Convert to a data frame so to add clinical columns

gene_expres_matrix_df <- as.data.frame(gene_expres_matrix)

# 5.4 Keep only ER+ patients

metadata_gse_2034_er_pos <- metadata_gse_2034  %>% 
  filter(ER_STAT == "ER+")

# 5.5 Keep oly patients that are found on the metadata

gene_expres_matrix_gse2034 <- gene_expres_matrix[,colnames(gene_expres_matrix) %in% metadata_gse_2034_er_pos$file_name] %>% 
  as.data.frame() %>% 
  t() 


saveRDS(gene_expres_matrix_gse2034, "./output_data/gene_expres_matrix_gse2034.RDS")
saveRDS(metadata_gse_2034_er_pos, "./output_data/metadata_gse_2034_er_pos.RDS")


rm(list = ls())
gc()
