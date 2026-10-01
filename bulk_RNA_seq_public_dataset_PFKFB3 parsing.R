###########################################################################################################################
## Project: Short term PFKFB3 mouse project
## code: public_dataset_1
## Author: Monojit
## Date: 20/02/2025
## Aim: Parsing out the PFKFB3 fold change in several inflammatory mediators
##########################################################################################################################

# Laod the libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(WGCNA)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DT)
library(data.table)
library(stats)
library(RColorBrewer) # for a colourful plot
library(CorLevelPlot)
library(ggrepel)


##################################################################################################
# Functions that are to be used in this code
# Id conversion from gene id to gene symbol
id_to_gene <- function(deg) {
  # Map ENTREZ IDs to gene symbols
  gene_syn <- mapIds(org.Hs.eg.db,
                     keys = deg$gene_id,
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multiVals = "first")
  
  # Add gene symbols to dataframe
  deg$gene_name <- gene_syn[match(deg$gene_id, names(gene_syn))]
  
  return(deg)
}
##################################################################################################

### For GSE118815
counts_data_a<- read.table("D:/Kroon_Lab/Analysis/Bulk_Cell/Public_dataset/GSE184512/GSE184512_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)

# get metadata --------
gse <- getGEO(GEO = 'GSE184512', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

metadata.modified <- metadata %>%
  dplyr::select(11,12) %>%
  setNames(c("Genotype", "Treatment")) %>%
  mutate(genotype = gsub("Genotype/variation:", "", Genotype)) %>%
  mutate(treatment = gsub("Treatment: ", "", Treatment)) %>%  
  mutate(group = gl(4, 2, labels = c("A","B", "C", "D")))



all(colnames(counts_data_a) %in% rownames(metadata.modified))
all(colnames(counts_data_a) == rownames(metadata.modified))


dds <- DESeqDataSetFromMatrix(countData = counts_data_a,
                              colData = metadata.modified,
                              design = ~ group)


dds_result <- DESeq(dds)

res_BA <- as.data.frame(results(dds_result, contrast = c("group", "B", "A"))) %>%
  mutate(
    gene_id = rownames(.)
    )

res_CA <- as.data.frame(results(dds_result, contrast = c("group", "C", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_DA <- as.data.frame(results(dds_result, contrast = c("group", "D", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_DB <- as.data.frame(results(dds_result, contrast = c("group", "D", "B"))) %>%
  mutate(
    gene_id = rownames(.))

########################################################################################################################

### For GSE131590
counts_data_a<- read.table("D:/Kroon_Lab/Analysis/Bulk_Cell/Public_dataset/GSE131590/GSE131590_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)

# get metadata --------
gse <- getGEO(GEO = 'GSE131590', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

metadata.modified <- metadata %>%
  dplyr::select(1, 8) %>%
  setNames(c("sample", "Treatment_Type"))

metadata.modified_m <- metadata.modified %>%
  mutate(Group = case_when(
    Treatment_Type == "LPS_sup" ~ "A",
    Treatment_Type == "Can_sup" ~ "B",
    Treatment_Type == "RPMI_sup" ~ "C",
    Treatment_Type == "RPMI_sup+IL1RA" ~ "D",
    Treatment_Type == "Can_sup+ IL1RA+ TNF_Ab" ~ "E",
    Treatment_Type == "LPS_sup+ IL1RA+ TNF_Ab" ~ "F",
    Treatment_Type == "Aspergillus fumigatus_24h" ~ "G",
    Treatment_Type == "Aspergillus fumigatus_4h" ~ "H",
    Treatment_Type == "Candida" ~ "F12",
    Treatment_Type == "Candida albicans_24h" ~ "J",
    Treatment_Type == "LPS" ~ "K",
    Treatment_Type == "Mycobacterium tuberculosis_4h" ~ "L",
    Treatment_Type == "Pseudomonas aeruginosa_24h" ~ "M",
    Treatment_Type == "RPMI_24h" ~ "N",
    Treatment_Type == "RPMI_4h" ~ "O",
    Treatment_Type == "Spneu" ~ "P",
    Treatment_Type == "Strep_sup" ~ "Q",
    Treatment_Type == "Streptococcus pneumoniae_4h" ~ "R",
    Treatment_Type == "TNF" ~ "S",
    Treatment_Type == "Candida albicans_4h" ~ "T",
    Treatment_Type == "IL-1b" ~ "U",
    Treatment_Type == "Mycobacterium tuberculosis_24h" ~ "V",
    Treatment_Type == "Pseudomonas aeruginosa_4h" ~ "W",
    Treatment_Type == "RPMI" ~ "Z",
    Treatment_Type == "Streptococcus pneumoniae_24h" ~ "A11",
    Treatment_Type == "Strep_sup+ IL1RA+ TNF_Ab" ~ "A12",
    TRUE ~ "Other"
  ))




all(colnames(counts_data_a) %in% rownames(metadata.modified_m))
all(colnames(counts_data_a) == rownames(metadata.modified_m))


dds <- DESeqDataSetFromMatrix(countData = counts_data_a,
                              colData = metadata.modified_m,
                              design = ~Group)


dds_result <- DESeq(dds)

res_S_Z <- as.data.frame(results(dds_result, contrast = c("Group", "S", "Z"))) %>%
  mutate(
    gene_id = rownames(.))

res_U_Z <- as.data.frame(results(dds_result, contrast = c("Group", "U", "Z"))) %>%
  mutate(
    gene_id = rownames(.))

res_K_Z <- as.data.frame(results(dds_result, contrast = c("Group", "K", "Z"))) %>%
  mutate(
    gene_id = rownames(.))


###############################################################################################################
### For GSE213111
counts_data_a<- read.table("D:/Kroon_Lab/Analysis/Bulk_Cell/Public_dataset/GSE213111/GSE213111_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)

# get metadata --------
gse <- getGEO(GEO = 'GSE213111', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)


ka <- c(rep(LETTERS[1], 6), rep(LETTERS[3:26], each=3))

metadata.modified <- metadata %>%
  dplyr::select(1, 13) %>%
  setNames(c("sample", "treatment")) %>%
  mutate(
    treatment = gsub("treatment: ", "", treatment)
  ) %>% mutate(group= ka)



all(colnames(counts_data_a) %in% rownames(metadata.modified))
all(colnames(counts_data_a) == rownames(metadata.modified))


dds <- DESeqDataSetFromMatrix(countData = counts_data_a,
                              colData = metadata.modified,
                              design = ~ group)


dds_result <- DESeq(dds)


res_X_A <- as.data.frame(results(dds_result, contrast = c("group", "X", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_M_A <- as.data.frame(results(dds_result, contrast = c("group", "M", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_Y_A <- as.data.frame(results(dds_result, contrast = c("group", "Y", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_P_A <- as.data.frame(results(dds_result, contrast = c("group", "P", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_Z_A <- as.data.frame(results(dds_result, contrast = c("group", "Z", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_L_A <- as.data.frame(results(dds_result, contrast = c("group", "L", "A"))) %>%
  mutate(
    gene_id = rownames(.))


###########################################################################################################
### For GSE118815
counts_data<- read.table("D:/Kroon_Lab/Analysis/Bulk_Cell/Public_dataset/GSE118815/GSE118815_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)

# get metadata --------
gse <- getGEO(GEO = 'GSE118815', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

metadata.modified <- metadata %>%
  as_tibble() %>%
  dplyr::select(2,10,11) %>%
  setNames(c("sample", "endmttreatment", "lvtreatment")) %>%
  dplyr::mutate(group = gl(6, 3, labels = c("A","B", "C", "D", "E", "F"))) %>%
  as.data.frame()

rownames(metadata.modified) <- metadata.modified$sample


all(colnames(counts_data) %in% rownames(metadata.modified))
all(colnames(counts_data) == rownames(metadata.modified))

#Differential expression
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metadata.modified,
                              design = ~ group)


dds_result <- DESeq(dds)

res_B_A <- as.data.frame(results(dds_result, contrast = c("group", "B", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_C_A <- as.data.frame(results(dds_result, contrast = c("group", "C", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_D_A <- as.data.frame(results(dds_result, contrast = c("group", "D", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_E_A <- as.data.frame(results(dds_result, contrast = c("group", "E", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_F_A <- as.data.frame(results(dds_result, contrast = c("group", "F", "A"))) %>%
  mutate(
    gene_id = rownames(.))

res_F_E <- as.data.frame(results(dds_result, contrast = c("group", "F", "E"))) %>%
  mutate(
    gene_id = rownames(.))

#######################################################################################################
# Parsing out the PFKFB3 gene Fold change

spc_gene_path<-function(rest,expset){
  # Map ENTREZ IDs to gene symbols
  gen_in_pway <- "PFKFB3"
  
  if (expset=="Public"){
    gene_syn <- mapIds(org.Hs.eg.db,
                       keys = rest$gene_id,
                       column = "SYMBOL",
                       keytype = "ENTREZID")
    # Add gene symbols to dataframe
    rest$gene_name <- gene_syn[match(rest$gene_id, names(gene_syn))]}
  
  if (expset == "Public"){
    rest_pased<-rest %>%
      filter(pvalue <.05 & padj <.1) %>%
      filter(gene_name %in% gen_in_pway) %>%
      dplyr::select(8, 2)
  } else if (expset=="Private") {
    rest_pased<-rest %>%
      filter(pvalue <.05 & padj <.1) %>%
      filter(gene_name %in% gen_in_pway) %>%
      select(1, 3)
  }
  
  return(rest_pased)
  
}

public_conditions <- list(
  "10 ng/ml TGFB2 + 1ng/ml IL1B (7days) vs unstim (GSE118815)"= res_D_A, 
  "25ng/mL TNF(12 h) vs Unstim (GSE184512)"= res_CA, 
  "10ng/ml TNF (6h) vs Unstim (GSE131590)"=res_S_Z, 
  "1000ng/ml LPS (6h) vs Unstim (GSE131590)"=res_K_Z, 
  "10ng/ml IL1B (6h) vs Unstim (GSE131590) "=res_U_Z,
  "10ng/ml TNF (24h) vs Unstim (GSE23111)"=res_X_A, 
  "10ng/ml IFN (24h) vs Unstim (GSE23111)"=res_Y_A, 
  "(each 10ng/ml)TNF + IFN (24h) vs Unstim (GSE23111)"=res_Z_A)


results_list_public<-lapply(public_conditions, function(x) spc_gene_path(x, expset="Public"))

results_df_public <- bind_rows(results_list_public, .id = "Conditions")

results_df_public <- results_df_public %>%
  mutate(FoldChange = ifelse(log2FoldChange > 0,
                             2^log2FoldChange,
                             -1/(2^abs(log2FoldChange)))) %>%
  
  arrange(FoldChange) %>%
  mutate(Conditions = factor(Conditions, levels = Conditions))

# visualization of PFKFB3 Fold change
ggplot(results_df_public,
       aes(x = Conditions,
           y = FoldChange)) +
  
  geom_bar(stat = "identity",
           fill = "grey",
           width = 0.65,
           color = "black",
           linewidth = 0.5) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = NULL,
    y = "Fold Change"
  ) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust=1.5, face = "bold"),
    axis.line = element_line(linewidth = 0.6),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  coord_flip()
