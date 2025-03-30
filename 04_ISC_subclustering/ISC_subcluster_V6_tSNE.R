#' ---
#' title: "ISC subclutering"
#' author: "Jingwei Song"
#' date: "March 30, 2025"
#' ---

setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/i-cell_subclustering")

#' load packages
#+ message = FALSE
library(Seurat)
library(cowplot) 
library(ggplot2)
library(dplyr)
load("../CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_final.RData")
DimPlot(seurat_integrated, label = T, label.size = 6)
# UMAP has two res, current res = 0.3

# How many cells now we have in C3 cluster?
colnames(seurat_integrated@meta.data)
table(seurat_integrated@meta.data$integrated_snn_res.0.3) #3457 cells
table(seurat_integrated@meta.data$integrated_snn_res.0.5) #3394
# lower resolution has more cells

# Subclustering

# extract just cluster 3
C3 <- subset(seurat_integrated, idents = "3")
C3 

DefaultAssay(C3) <- "RNA"

#' integrate CCA 
#+ message = FALSE  
C3.list <- SplitObject(C3, split.by = "sample")
C3.list # Live_cells, Methanol_New

#' perform SCTransform individually
#+ message = FALSE, warning = FALSE
C3.list <- lapply(C3.list, FUN = function(x){
  x <- SCTransform(x, vst.flavor = "v2") # v2 instead of default
})

# 3000 maybe not necessary
integ_features <- SelectIntegrationFeatures(object.list = C3.list, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
C3.list <- PrepSCTIntegration(object.list = C3.list, 
                              anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = C3.list, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

seurat_integrated #3457 cells 

#' Run PCA
#+ message = FALSE
seurat_integrated <- RunPCA(seurat_integrated, npcs = 25)

seurat_t <- RunTSNE(seurat_integrated)

seurat_t <- FindNeighbors(seurat_t, dims = 1:25) %>%
  FindClusters(resolution = 0.2)
TSNEPlot(seurat_t, label = T, label.size = 6,
         split.by = "sample")
TSNEPlot(seurat_t, label = T, label.size = 6)


# res from 0.3 to 0.2, 7 cluster become 6 clusters.

DefaultAssay(seurat_t) <- "RNA"
#' C0, Piwi, NOP56,NOP58,GNL3
FeaturePlot(seurat_t, slot = "data", features = c("HyS0050.7","HyS0073.68", "HyS0155.10","HyS0059.86"), reduction = "tsne",
            order = T, max.cutoff = 'q95')

#'  C1, HyS0022.107 (NR2E1), HyS0008.263 (Ncol1), PCNA
FeaturePlot(seurat_t, slot = "data", features = c("HyS0022.107","HyS0008.263", "HyS0061.60"), reduction = "tsne",
            order = T, max.cutoff = "q95")

#' C2
# HyS0034.90, Hippocalcin-like; HyS0085.53, ELAV, HyS0012.308 Sox22, 
FeaturePlot(seurat_t, slot = "data", features = c("HyS0034.90","HyS0085.53","HyS0012.308"), reduction = "tsne",
            order = T, max.cutoff = "q95")

#' C3
# Fat4, Chit1, astacin-2, delta2 (notch signaling),notch-b (HyS0001.72)

FeaturePlot(seurat_t, features = c("HyS0030.178","HyS0041.99","HyS0013.50", "HyS0088.7"), reduction = "tsne",
            order = T, max.cutoff = "q95")


#' C4
#HyS0013.62, HyTSR1, Myc
FeaturePlot(seurat_t, slot = "data", features = c("HyS0013.62","HyS0005.84"),
            order = T,  max.cutoff = "q95")


#' where are these subcluster cells on the large UMAP
#' 
load("../CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_final.RData")
i <- NA
for (i in 0:5) {
  Iwant <- WhichCells(seurat_t, idents = i)
  print(DimPlot(seurat_integrated, cells.highlight = Iwant) + labs(title = paste0('Cluster', i)))
}

# C5 is likely an artefact, not in the C3 of large UMAP. 
# C3 subcluster is likely real, but what are they? 
seurat_t <- subset(seurat_t, idents = "5", invert = T)
DimPlot(seurat_t, label = T)
#save(seurat_t, file = "ISC_subcluter_Seurat_V6_tnse.RData")

#' plot multiple genes at once
#' Load Seurat and patchwork
library(Seurat)
library(patchwork)

#' List of six genes to plot
genes_to_plot <- c("HyS0050.7", "HyS0073.68", "HyS0022.107", "HyS0061.60", "HyS0028.205", "HyS0012.308")

#' Generate FeaturePlots for each gene
plot_list <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(seurat_t, features = gene, order = T, max.cutoff = 'q95')
})

# Combine the plots into a 2x3 grid using patchwork
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 3)

# Display the combined plot
print(combined_plot)


sessionInfo()


