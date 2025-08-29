#load required packages
library(Seurat)
library(tidyverse)
library(BiocManager)
library(glmGamPoi)

#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")

# see UMAP
DimPlot(seurat_integrated, label = F, label.size = 5, pt.size=1)

# Subset clusters 7 and 14
sub_cluster7_14_resolution0.1 <- subset(seurat_integrated, idents = c(7, 14))

# Save original cluster IDs
sub_cluster7_14_resolution0.1$orig_cluster <- Idents(sub_cluster7_14_resolution0.1)

# Switch back to SCT assay
DefaultAssay(sub_cluster7_14_resolution0.1) <- "SCT"

# Re-run SCT on the subset
sub_cluster7_14_resolution0.1 <- SCTransform(sub_cluster7_14_resolution0.1, verbose = FALSE,
                                             vst.flavor = "v2")

# Dimensional reduction and clustering
sub_cluster7_14_resolution0.1 <- RunPCA(sub_cluster7_14_resolution0.1)
sub_cluster7_14_resolution0.1 <- FindNeighbors(sub_cluster7_14_resolution0.1, dims = 1:15)
sub_cluster7_14_resolution0.1 <- FindClusters(sub_cluster7_14_resolution0.1, resolution = 0.1)
sub_cluster7_14_resolution0.1 <- RunUMAP(sub_cluster7_14_resolution0.1, dims = 1:15)

# Visualize
DimPlot(sub_cluster7_14_resolution0.1, label = TRUE, pt.size = 1)


# Visualize ORIGINAL clusters (7 vs 14)
DimPlot(sub_cluster7_14_resolution0.1, group.by = "orig_cluster", pt.size = 1) + 
  ggtitle("Original clusters (7 vs 14)")


######
#Getting DE lists for my clusters

DefaultAssay(sub_cluster7_14_resolution0.1) <- "RNA"

# Differential expression for all reclustered groups
de_markers <- FindAllMarkers(
  sub_cluster7_14_resolution0.1,
  only.pos = TRUE,          # only positive markers
  min.pct = 0.1,           # expressed in at least 25% of cells in a cluster
  min.diff.pct = 0.5,
  logfc.threshold = 1    # minimum log2 fold-change
)

#read in gene annotation files
annotation <- read.csv("~/Desktop/R files/hsym_combined_annotations_2.csv")
head(annotation)

# Need to reformat first column
library(stringr)

Iwant <- str_replace(annotation$gene_id, "Hsym[|]", "")
head(Iwant)
annotation$gene_id <- Iwant
head(annotation)

#de.markers <- rownames_to_column(de.markers, var = "gene_id")
head(de_markers)
head(annotation)

# tmp <- left_join(de.markers, annotation, by = "gene_id") # de.markers column 7 change name
colnames(de_markers)[7] <- "gene_id"
colnames(de_markers)

# Now it works
tmp <- left_join(de_markers, annotation, by = "gene_id")

head(tmp)
#View(tmp)


write.csv(tmp, "subcluster7_14_DE_markers_resolution01.csv", row.names = TRUE)


