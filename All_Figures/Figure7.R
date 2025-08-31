#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final.Rdata")
# see UMAP
DimPlot(seurat_integrated, label = F, label.size = 5, pt.size=1)
# need to run this to change the default assay for the current UMAP
DefaultAssay(seurat_integrated) <- "SCT"


tmp <- seurat_integrated
table(Idents(tmp))
levels(x = tmp)

cluster.ids <- c("epithelial cells", "nematocytes", "nematocytes", "i-cells", "nematocytes",
                 "epithelial cells", "nematoblasts", "neurons", "nematocytes", "immune cells",
                 "gland cells", "gland cells", "gland cells", "gland cells", "neurons",
                 "nematoblasts", "nematoblasts", "epithelial cells", "epithelial cells") 
names(cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, cluster.ids)

a <- c('lightgrey', 'lightgrey', '#0247FE', 'lightgrey', 'lightgrey','lightgrey', 'lightgrey')

DimPlot(tmp, reduction = "umap", label.size = 5, pt.size = 1.0, cols = a) + 
  NoLegend() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "cm")
  )

library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots

#load Jing-wei's new i-cell subcluster saved on desktop_v6
#This is now the version without c5
load("~/Desktop/R files/ISC_subcluter_Seurat_V6_tnse.RData")
# see UMAP
DimPlot(seurat_t, label = FALSE, label.size = 5, pt.size = 2) + 
  labs(title = NULL) +  
  theme(
    axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
    axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
    axis.title.y = element_text(face = "bold", size = 30),
    axis.title.x = element_text(face = "bold", size = 30),
    axis.line = element_line(linewidth = 2.5),
    axis.ticks = element_line(linewidth = 2.5),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = "none"  # This removes the legend
  )

#HyS0050.7_Piwi1
FeaturePlot(seurat_t, slot = "data", "HyS0050.7", reduction = "tsne",
            order = T, max.cutoff ="q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

#HyS0073.68_NOP56
FeaturePlot(seurat_t, slot = "data", "HyS0073.68", reduction = "tsne",
            order = T, max.cutoff = "q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0022.107_NR2E1
FeaturePlot(seurat_t, slot = "data", "HyS0022.107", reduction = "tsne",
            order = T, max.cutoff = "q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0061.60_PCNA
FeaturePlot(seurat_t, slot = "data", "HyS0061.60", reduction = "tsne",
            order = T, max.cutoff = "q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#HyS0028.205_Neurogenin
FeaturePlot(seurat_t, slot = "data", "HyS0028.205", reduction = "tsne",
            order = T, max.cutoff = "q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#HyS0012.308_Sox22
FeaturePlot(seurat_t, slot = "data", "HyS0012.308", reduction = "tsne",
            order = T, max.cutoff = "q95", pt.size =3, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 
