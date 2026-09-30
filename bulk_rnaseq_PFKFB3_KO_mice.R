######################################################################################################################
### Bulk RNA-seq of endothelial cells isolated from PFKFB3 KO and control mice
## Author: Monojit
## Date: 18/10/2025
## samples: PFKFB3 Knockout vs Control mice
#######################################################################################################################

# Set working path
setwd("D:/Kroon_Lab/Analysis/Bulk_Cell/PFKFB3")


# Load required packages
library(readr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)


# Read the data
data <- read.csv("gene_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

dim(data)
head(data)

# Extract count matrix and gene lengths
rownames(data) <- data$Ensembl_ID

# Extract gene lengths and convert to kilobases
gene_lengths <- data$Length
names(gene_lengths) <- data$Ensembl_ID
gene_lengths_kb <- gene_lengths / 1000  # Convert to kilobases

# TPM normalization function
TPM_normalization <- function(counts_matrix, gene_lengths_kb) {
  rpk <- sweep(counts_matrix, 1, gene_lengths_kb, FUN = "/")
  scaling_factors <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")
  return(tpm)
}

# Run TPM normalization
tpm_matrix <- TPM_normalization(data[,3:4], gene_lengths_kb)
tpm_df <- as.data.frame(tpm_matrix)
tpm_df$gene_id <- sub("\\..*", "", rownames(tpm_df))
head(tpm_df)


# Map the ENSEMBL -> SYMBOL
gene_name <- mapIds( org.Mm.eg.db, 
                     keys = tpm_df$gene_id, 
                     column = "SYMBOL",
                     keytype = "ENSEMBL", 
                     multiVals = "first")

tpm_df$gene_name <- gene_name[match(tpm_df$gene_id, names(gene_name))]

# Removing the NA and taking mean of same gene
tpm_df_agg <- tpm_df %>%
  na.omit() %>%
  group_by(gene_name) %>%
  summarise(
    PFKFB3_WT = mean(PFKFB3_WT, na.rm = TRUE),
    PFKFB3_KO = mean(PFKFB3_KO, na.rm = TRUE),
  ) %>%
  as.data.frame()


tpm_df_agg$diff <- tpm_df_agg$PFKFB3_KO - tpm_df_agg$PFKFB3_WT
rownames(tpm_df_agg) <- toupper(tpm_df_agg$gene_name)

# heatmap of cam genes
cam <- c( "ICAM2", "CDH2", "CDH5", "CDH1",
          "PECAM1", "ICAM1", "ICAM3", "VCAM1", "NCAM1",
          "ITGAL", "ITGAM", "ITGAX", "ITGA4", "ITGB1", "ITGB2", "ITGB7", 
          "SELE", "SELL", "SELP", 
          "F11R", "JAM2", "JAM3", "SLC35C1")

tpm_cam<-tpm_df_agg[cam,] %>%
  as.data.frame() %>%
  na.omit

rownames(tpm_cam) <- tpm_cam$gene_name

log_tpm_cam <- log2(tpm_cam[,2:3] + 1)
colnames(log_tpm_cam) <- c("EC-PFKFB3-WT", "EC-PFKFB3-KO")

p2 <- pheatmap::pheatmap(
  t(log_tpm_cam),
  scale = "column",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  col = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(-1, 1, length.out = 100),
  border_color = NA,
  treeheight_col = 0,
  treeheight_row = 0,
  fontsize_row = 12,
  fontsize_col = 12,
  cellheight = 80,
  main = "The Log2(TPM + 1) Expression of Leukopcyte Adhesion Genes",
  angle_col = 45
)


# Add padding between title and heatmap
p2$gtable$heights[1] <- unit(1.5, "cm")

final_cam<- gridExtra::grid.arrange(
  p2$gtable,
  vp = grid::viewport(width = 0.9, height = 1)
)


# Heatmap expression of endomt genes

endomt<- c(
  "NOTCH1", "SNAI1", "KLF4", "LEF1", "TCF7", "DLL4", "JAG1", "TNF", "IL1B", "SMAD5", "CDH2", "CDH5", "TEK",
  "ACTA1", "COL1A1", "COL3A1", "FN1", "VIM", "WNT3A", "WNT5A", "S100A4", "TGFB1", "TGFB2", "TGFBR1", "SMAD2", "SMAD3", "SMAD7",
  "ZEB1", "ZEB2", "SNAI2", "TWIST1", "HES1", "EZH2", "HDAC3", "BRD4", "NOTCH4", "PECAM1", "VWF", "CLDN5", "TIE1", 
  "RELA", "NFKB1", "HIF1A"
  
)

tpm_endo<-tpm_df_agg[endomt,] %>%
  as.data.frame() %>%
  na.omit

rownames(tpm_endo) <- tpm_endo$gene_name

log_tpm_endo <- log2(tpm_endo[,2:3] + 1)
colnames(log_tpm_endo) <- c("EC-PFKFB3-WT", "EC-PFKFB3-KO")

p4 <- pheatmap::pheatmap(
  t(log_tpm_endo),
  scale = "column",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  col = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(-1, 1, length.out = 100),
  border_color = NA,
  treeheight_col = 0,
  treeheight_row = 0,
  fontsize_row = 12,
  fontsize_col = 12,
  cellheight = 80,
  main = "The Log2(TPM + 1) Expression of EndoMT Genes",
  angle_col = 45
)


# Add padding between title and heatmap
p4$gtable$heights[1] <- unit(1.5, "cm")

final_endo<- gridExtra::grid.arrange(
  p4$gtable,
  vp = grid::viewport(width = 0.9, height = 1)
)

############################################################################################################################################
#Figures in the paper

## Figure 8H:
final_cam

## Figure 4A
final_endo


