## Nov 4, 2022 ##
# read in raw matrices for all 8 ACME samples
# updated Mar 2, 2023 by jingwei

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/10X/Jan_2023")

# Create a Seurat object for each sample
for (file in c("BR4_lib1_filtered_feature_bc_matrix", 
               "BR4_lib2_filtered_feature_bc_matrix",
               "BR5_lib1_filtered_feature_bc_matrix",
               "BR5_lib2_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Merge all samples into a single Seurat object
merged_seurat <- merge(x = BR4_lib1_filtered_feature_bc_matrix, 
                       y = c(BR4_lib2_filtered_feature_bc_matrix,
                             BR5_lib1_filtered_feature_bc_matrix,
                             BR5_lib2_filtered_feature_bc_matrix),
                       add.cell.id = c("BR4_lib1", "BR4_lib2", "BR5_lib1", "BR5_lib2"))

# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

### QC ###

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^HyS0613.")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Create sample column
metadata$sample <- NA
metadata$sample <- rep("Methanol_New", dim(metadata)[1])

# Rename columns, "new.name = old.name"
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident, # seq_folder 
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

tail(merged_seurat@meta.data)

# OPTIONAL:
#save(merged_seurat, file="data/merged_seurat_step01.RData")

# Visualize the number of cell counts per sample
merged_seurat@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via histogram
merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
merged_seurat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
merged_seurat@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 300 ) & 
                            (nGene >= 100 ) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Gene level filtering
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

### 2/13/23 #####
# run until line 71

# This cutoff to match with ACME
filtered_seurat <- subset(x = filtered_seurat, 
                          subset= (nUMI >= 300 & nUMI <=100000) & 
                            (nGene >= 100 & nGene <7000) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

filtered_seurat <- SCTransform(filtered_seurat, vst.flavor = "v2") %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution = 0.5)  

DimPlot(filtered_seurat,label = T) + NoLegend() + NoAxes() 
# similar to what's seen before
# suspect C3 is artefact cluster

FeaturePlot(filtered_seurat, features = c("HyS0018.100", "HyS0061.60", #HMG, PCNA
                                          "HyS0008.263", "HyS0030.203", #Ncol1, Nematocilin
                                          "HyS0041.99", "HyS0004.446", #chitinase,mucs
                                          "HyS0085.53", "HyS0013.338", #ELAV, Rfamide
                                          "HyS0017.253", "HyS0010.323",#Wnt5, Pitx
                                          "HyS0001.667", "HyS0011.54"), 
            order = T, min.cutoff = 'q10')

#save(filtered_seurat, file="output/filtered_seurat_PC50_dims15_res0.5_w_artefact.RData")

# removed C3
no_C3 <- subset(filtered_seurat, idents = "3", invert = T)

DimPlot(no_C3)

# rerun UMAP and recluster
no_C3 <- SCTransform(no_C3, vst.flavor = "v2") %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution = 0.5)  

DimPlot(no_C3,label = T) + NoLegend() + NoAxes()

FeaturePlot(no_C3, features = c("HyS0018.100", "HyS0061.60", #HMG, PCNA
                                          "HyS0008.263", "HyS0030.203", #Ncol1, Nematocilin
                                          "HyS0041.99", "HyS0004.446", #chitinase,mucs
                                          "HyS0085.53", "HyS0013.338", #ELAV, Rfamide
                                          "HyS0017.253", "HyS0010.323",#Wnt5, Pitx
                                          "HyS0001.667", "HyS0011.54"), 
            order = T, min.cutoff = 'q10')

#save(no_C3, file="output/filtered_seurat_PC50_dims15_res0.5_wo_artefact.RData")