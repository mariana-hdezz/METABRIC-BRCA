library(clusterProfiler)
library(org.Hs.eg.db)

proof_genes_surv <- readRDS("./output_data/proof_genes_surv.RDS")
proof_genes_rec <- readRDS("./output_data/proof_genes_rec.RDS")

for (i in c(1, 2)) {
  if(i == 1){
    analysis <- "recurrence"
    proof_genes <- proof_genes_rec
  }else{
    analysis <- "survival"
    proof_genes <-proof_genes_surv
  }


# 1. Your list of genes
genes_to_test <- proof_genes
genes_to_test <- toupper(trimws(genes_to_test))

# Convert Symbols to Entrez IDs
gene_conv <- bitr(genes_to_test, 
                  fromType = "SYMBOL", 
                  toType   = "ENTREZID", 
                  OrgDb    = org.Hs.eg.db)

# Run enrichment using the Entrez IDs
go_results <- enrichGO(gene = gene_conv$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = 'ENTREZID', 
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       readable = TRUE) #

if(i == 1){
  go_results_rec <- go_results
}else{
  go_results_surv <- go_results
}

}

as.data.frame(go_results_rec)
