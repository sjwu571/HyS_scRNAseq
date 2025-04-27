# Differential Expression and Cluster Annotations 
# updated Apr 23, 2025 by Jingwei

setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/DE")

library(Seurat)
library(cowplot) 
library(ggplot2)
library(dplyr)
load("../Seurat_integrated_features3k_PC50_dims20_res0.5_final.RData")
DimPlot(seurat_integrated, label = T, label.size = 6)

####### DE ###########

DefaultAssay(seurat_integrated) <- "RNA" # used integrated before

de.markers <- FindAllMarkers(seurat_integrated, only.pos = T, 
                             min.pct = 0.3, logfc.threshold = 1)

#read in gene annotation files
annotation <- read.csv("../../i-cell_subclustering/DE/hsym_combined_annotations_2.csv")
head(annotation)

# Need to reformat first column
library(stringr)

Iwant <- str_replace(annotation$gene_id, "Hsym[|]", "")
head(Iwant)
annotation$gene_id <- Iwant
head(annotation)

#de.markers <- rownames_to_column(de.markers, var = "gene_id")
head(de.markers)
head(annotation)

# tmp <- left_join(de.markers, annotation, by = "gene_id") # de.markers column 7 change name
colnames(de.markers)[7] <- "gene_id"
colnames(de.markers)

# Now it works
tmp <- left_join(de.markers, annotation, by = "gene_id")

head(tmp)
#View(tmp)

# write to home directory
write.csv(tmp, file = "DE_markers_RNA_minpct0.3.csv")

# Check some top DE genes for each subcluster

# sessionInfo()
