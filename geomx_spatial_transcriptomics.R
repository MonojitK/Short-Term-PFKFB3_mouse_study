######################################################################################################################
### Spatial Transcriptomics ( GeoMX ) of Carotid Plaques
## Author: Monojit
## Date: 18/08/2025
######################################################################################################################
## Samples: PFKFB3 EC specific KO (PFKFB3fl/flCdh5Cre+/-LDLR-/- ) vs control (PFKFB3fl/flCdh5Cre-/-LDLR-/- )
## Sample info: 
##             PFKFB3 EC specific KO ( Mouse ID ): A42, A50, A52, A57
##             control ( Mouse ID ): A33, A34, A39, A43
#######################################################################################################################
## Dataset: Fastqc is converted to the resulting datasets using GeoMX NGS pipeline by nanostring
## Resulting datasets: 
##                  1. DCC file
##                  2. pKC file from Nanonstring ( Mouse whole Transcriptome Atlas)
##                  3. Annotation file from Lab worksheet ( Need to modify the Lab work sheet)
## Resulting datasets:
##                  1. Counts data
##                  2. metadata
##                  3. Gene metadata
#########################################################################################################################
# =======================================
# 2. Define file paths
# =======================================
# Working Directory
setwd("D:/Kroon_Lab/Analysis/Spatial_Transcriptomics/GeoMX/Analysis")

# Data Directory
datadir <- file.path("D:/Kroon Lab/Spatial Transcriptomics/Spatial Transcriptomics/Analysis/Result_GeoMX_NGS_pipeline")

# Laoding the libraries
# =======================================
# 1. Load libraries
# =======================================
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(biomaRt)
library(readr)
library(tidyverse)
library(dplyr)
library(SpatialExperiment)
library(standR)
library(ggplot2)
library(ggalluvial)
library(pheatmap)

DCCFiles <- dir(
  path = file.path("Result_GeoMX_NGS_pipeline", "Final_Counts", "DCC-20250818"),
  pattern = ".dcc$",
  full.names = TRUE,
  recursive = TRUE
)

PKCFiles <- file.path("Result_GeoMX_NGS_pipeline","Input_raw_metadata","Mm_R_NGS_WTA_v1.0.pkc")
SampleAnnotationFile <- file.path("Result_GeoMX_NGS_pipeline","Input_raw_metadata","annotations.xlsx")

# =======================================
# 3. Load GeoMx DCC data
# =======================================
chr4data <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = "250428_2302_AR_spatial_20250429",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("Aoi", "ROI_whole"),
  experimentDataColNames = c("Panel")
)

# =======================================
# 4. Extract counts and annotations
# =======================================
countFile <- as.data.frame(chr4data@assayData$exprs)
sampleAnnoFile <- as.data.frame(chr4data@phenoData@data)
featureAnnoFile <- as.data.frame(chr4data@featureData@data)

# =======================================
# 5. Remove unwanted genes (optional)
# =======================================
# 1. Define genes and sample to remove
removeGenes <- c("D830030K20Rik", "Gm10406", "LOC118568634")

# 2. Create logical index of genes to keep
keepGenes <- !rownames(countFile) %in% removeGenes

# 3. Subset count matrix (remove unwanted genes and sample)
countFile <- countFile[keepGenes, ]

# 4. Subset feature annotation file (remove unwanted genes)
featureAnnoFile <- featureAnnoFile[keepGenes, ]

#sampleAnnoFile <- sampleAnnoFile[!(rownames(sampleAnnoFile) %in% "DSP-1001660020570-E-A01")]


featureAnnoFile$TargetName <- make.unique(as.character(featureAnnoFile$TargetName))

names_col <- sapply(strsplit(colnames(countFile), "-"), function(x) x[length(x)])
colnames(countFile) <- sub("\\.dcc$", "", names_col)
sampleAnnoFile$SegmentDisplayName <- colnames(countFile)
sampleAnnoFile$Slide_Name <- sampleAnnoFile$`slide name`
sampleAnnoFile <- sampleAnnoFile %>%
  mutate(Group = ifelse(Group == "Knock Out", "EC-PFKFB3-KO",
                        ifelse(Group == "Control", "EC-PFKFB3-WT", Group)))

countFile$TargetName <- featureAnnoFile$TargetName



# Adding gene length to feature annotation file
############################################################################################################################
# Connect to Ensembl mouse dataset
mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

# Get gene length info for your genes (TargetName)
gene_list <- featureAnnoFile$TargetName
gene_info <- getBM(
  attributes = c("external_gene_name", "gene_biotype", "transcript_length"),
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
)

gene_lengths <- gene_info %>%
  group_by(external_gene_name) %>%
  summarize(gene.lengths = max(transcript_length)) 


featureAnnoFile <- featureAnnoFile %>%
  left_join(gene_lengths, by=c("TargetName" = "external_gene_name"))


#write_tsv(counts,'counts.txt')
#write_tsv(metadata,'metadata.txt')
#write_tsv(genemeta,'genemeta.txt')

# Cellular deconvolution

library(SpatialDecon)
library(speckle)


spe<-readGeoMx(
  countFile[-1],
  sampleAnnoFile,
  featureAnnoFile,
  rmNegProbe = FALSE,
  colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
  coord.colnames = c("ROI_Coordinate_X", "ROI_Coordinate_Y")
)


# QC for sample level metadata
#########################################################################################################################
speq<-spe
colData(speq)$Slide_Name <- colData(speq)$'slide name'
colData(speq)$cell_id <- rownames(colData(speq))
colData(speq)$Group <- factor(
  colData(speq)$Group,
  levels = c("EC-PFKFB3-WT", "EC-PFKFB3-KO")
)

# Gene level QC
speq <- addPerROIQC(spe, sample_fraction = .9, min_count = 5, rm_genes = TRUE)
dim(speq)
metadata(speq) |> names()

supp_1<-plotGeneQC(speq, ordannots = "Slide_Name", col = Slide_Name, point_size = 2)

# ROI level QC
speq <- speq[, !colnames(speq) %in% "A01"]
colData(speq)$Group <- factor(colData(speq)$Group, levels = c("EC-PFKFB3-WT","EC-PFKFB3-KO"))
supp_2<-plotROIQC(speq, x_axis ="Nuclei", y_axis ="lib_size", x_threshold = 50, color = Group)

# Filter good quality ROI
qc <- colData(speq)$Nuclei > 50
table(qc)
spe <- speq[, qc & !is.na(qc)]

# Relative log expression distribution
a <- plotRLExpr(spe, ordannots = "Slide_Name", color = Slide_Name)
drawPCA(speq, assay = 2, color = Slide_Name)



#Dimensity Reduction
#PCA
set.seed(100)
speq <- scater::runPCA(speq)
pca_results <- reducedDim(speq, "PCA")
b<-drawPCA(speq, precomputed = pca_results, col = Slide_Name)
c<-drawPCA(speq, precomputed = pca_results, col = Roi_sample)
plotScreePCA(speq, precomputed = pca_results)
plotPairPCA(speq, col = Roi_sample, precomputed = pca_results, n_dimension = 4)

# Normalization 
spe_tmm <- geomxNorm(spe, method = "TMM", log = TRUE)
supp4 <- plotRLExpr(spe_tmm, assay = 2, color = Slide_Name) + ggtitle("TMM")

spe_tmm <- scater::runPCA(spe_tmm)
pca_results_tmm <- reducedDim(spe_tmm, "PCA")
g<-drawPCA(spe_tmm, precomputed = pca_results_tmm, col = Group)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Group)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Slide_Name)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = rownames(colData(spe_tmm)))

#QC statistics
neg_probe <- rownames(spe)[grepl("Neg", rownames(spe), ignore.case = TRUE)]
spd <- prepareSpatialDecon(spe_tmm,   
                           assay2use = "logcounts",
                           negProbeName = neg_probe, 
                           pool = rep(1, nrow(assay(spe_tmm, "logcounts"))))

rownames(spd$normCount) <- toupper(rownames(spd$normCount))
rownames(spd$backGround) <- toupper(rownames(spd$backGround))


# loading Plaque atlas signature
profile_matrix <- read.csv("../plaque_atlas_matrix.csv", row.names = 1)


rownames(profile_matrix) <- toupper(rownames(profile_matrix))
safeHT <- as.matrix(profile_matrix)

# Heatmap of the conserved markers from plaque atlas
supp_5<-heatmap(sweep(safeHT, 1, apply(safeHT, 1, max), "/"),
        labRow = NA, margins = c(10, 5))


res <- spatialdecon(norm = spd$normCount,
                    bg = spd$backGround,
                    X = safeHT,
                    align_genes = TRUE)

samples_subset <- colnames(spe_tmm)[colData(spe_tmm)$'slide name' %in%  c("33", "39_1", "34_1", "42", "43", "50", "52", "57")]

subset_prop <- res$prop_of_all[,samples_subset]

spe_sub <- spe_tmm[,samples_subset]
colnames(colData(spe_sub))[colnames(colData(spe_sub)) == "slide name"] <- "slide_name"

# Cell type decomposition based on slides
subset_prop %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  ggplot(aes(samples, prop, fill = CellTypes)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7)  +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top")+
  labs(fill = "Cell Types", y = "Proportion")

# finding the proportion of cells 
propslist <- convertDataToList(subset_prop, 
                               data.type = c("proportions"),
                               transform="asin",
                               scale.fac=colData(spe_sub)$Nuclei)

# Differential expression of the proportion of de-convoluted cell type
design <- model.matrix(~ 0 + Group, data = as.data.frame(colData(spe_sub)))

colnames(design) <- make.names(colnames(design))
colnames(design)

contr <- limma::makeContrasts(GroupEC.PFKFB3.KO  - GroupEC.PFKFB3.WT ,levels=design)


outs <- propeller.ttest(propslist, design, contr, robust=TRUE,trend=FALSE, sort=TRUE)

# Significant differential proportion
diff_ct <- outs %>% 
  filter(P.Value < 0.05) %>%
  rownames()

colData(spe_sub)$samples_id <- rownames(colData(spe_sub))
colData(spe_sub)$id <- rownames(colData(spe_sub))

all_cells <- subset_prop[rownames(outs),] %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  left_join(as.data.frame(colData(spe_sub)), by = c("samples"="samples_id")) 

all_cells$cell_count <- all_cells$prop * all_cells$Nuclei

#sub setting Macrophage, Fibroblast, Smooth Muscle Cell, EC differential proportion
spe_cells<- all_cells[all_cells$CellTypes %in% c("Macrophage", "Fibroblast", "Smooth.Muscle.Cell", "EC"), ]

# Violin Plot
pd <- position_dodge(width = 0.8)
ct2 <- spe_cells %>%
  ggplot(aes(x = Group, y = prop, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE, position = pd) +
  geom_boxplot(
    fill="white",
    width = 0.12,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black",
    position = pd
  ) +
  facet_wrap(~CellTypes, ncol = 4, nrow = 1) +
  scale_fill_manual(values = c(
    "EC-PFKFB3-KO" = "red",
    "EC-PFKFB3-WT" = "blue"
  )) +
  theme_bw(base_size = 14) +
  xlab("") +
  ylab("Proportion") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




# Stacked bar plot
ct1 <- spe_cells %>% 
  ggplot(aes(samples, prop, fill = CellTypes)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  geom_hline(yintercept = .96, linetype = "dashed", color = "red")+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  labs(fill = "Cell Types", y = "Proportion")


total_cell <- all_cells %>%
  group_by(Group, CellTypes) %>%
  summarise(
    mean_total_prop = mean(prop, na.rm = TRUE) * 100,     
    sd_total_prop   = sd(prop, na.rm = TRUE) * 100,       
    n               = n(),
    se_total_prop   = sd_total_prop / sqrt(n),            
    .groups = "drop"
  ) %>%
  mutate(
    mean_SE_label = paste0(round(mean_total_prop, 1), " ± ", round(se_total_prop, 1), "%")
  )



total_population <- all_cells %>%
  filter(CellTypes %in% c("Macrophage", "Smooth.Muscle.Cell", "EC", "Fibroblast")) %>%
  group_by(Group, samples) %>%
  summarise(total_prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  summarise(
    mean_total = mean(total_prop, na.rm = TRUE) * 100,  
    sd_total   = sd(total_prop, na.rm = TRUE) * 100,     
    n          = n(),
    se_total   = sd_total / sqrt(n),                     
    mean_SE_label = sprintf("%.1f ± %.1f%%", mean_total, se_total),
    .groups = "drop"
  )




# visualization used in the paper
###############################################################################
## Figure 3E
ct2

### Supplementary figure 1
supp_1

### Supplementary figure 2
supp_2

### Supplementary Figure 3
a+b+c

### Supplementary Figure 4
supp4+g


