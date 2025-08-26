##### sperm progenitors ######
library(Seurat)
library(ggplot2)
library(dplyr)
# find sperm progenitors in the live cell atlas
# load("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/Justin2019/aggr_fltd_PC50_dim15_res.0.5_v2.RData")
# DimPlot(aggr_fltd, label = T) # no, the sperm were filtered

# get sperm from Justin's UMAP
load("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Seurat/Justin's single cell/aggr_fltd.RData")

FeaturePlot(aggr_fltd, features = c("HyS0027.170","HyS0001.110")) 

# C6 is the sperm progenitor

sperm_prog <- subset(aggr_fltd, idents = "6")

sperm_prog # 394 cells

# need to run SCT before integration
sperm_prog <- SCTransform(sperm_prog, vst.flavor = 'v2')

# find C3 in the big UMAP
load("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")

ISC_prog <- subset(seurat_integrated, idents = "3")

ISC_prog # 3457 cells

# integrate both above
DefaultAssay(ISC_prog) <- "RNA"

ISC_prog <- SCTransform(ISC_prog, vst.flavor = 'v2')

split_seurat <- list(sperm_prog, ISC_prog)

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Run PCA (from ISC subclutering html)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 25)

seurat_t <- RunTSNE(seurat_integrated)
library(dplyr)
seurat_t <- FindNeighbors(seurat_t, dims = 1:25) %>%
  FindClusters(resolution = 0.2)

TSNEPlot(seurat_t, label = T)

# highlight cells expressing sperm marker
DefaultAssay(seurat_t) <- "RNA"

p1 <- FeaturePlot(seurat_t, features = c("HyS0027.170"), order = T)

# these sperm progenitors are in a specific cluster, C3

# highlight those cells in the ISC atlas
# get cell names for C3 subcluster in the ISC atlas
# first, give the tsne with sperm a new name
seurat_t_sperm <- seurat_t
# the other seurat_t (ISC UMAP in the paper)
load("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/i-cell_subclustering/ISC_subcluter_Seurat_V6_tnse.RData")
DimPlot(seurat_t, label = T)

Iwant <- WhichCells(seurat_t, idents = '3')

p2 <- TSNEPlot(seurat_t_sperm, cells.highlight = Iwant) + NoLegend()

# The unknown C3 in ISC atlas and sperm progenitors did not overlap.

p1+p2
