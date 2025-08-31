###Figure 2a panels

library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots

#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")
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
                 "nematoblasts", "nematoblasts", "epithelial cells", "epithelial cells", "nematocytes") 
names(cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, cluster.ids)

a <- c('lightgrey', '#66B032', '#0247FE', '#66B032', 'lightgrey','lightgrey', 'lightgrey')

DimPlot(tmp, reduction = "umap", order = TRUE, pt.size = 1, cols = a) + NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))




#load packages
library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots


#load Cnidogenesis UMAP saved on desktop
load("~/Desktop/R files/integrated_subset.RData")
# need to run this to change the default assay for the current UMAP
DefaultAssay(integrated_subset) <- "RNA"

#Fig2a panel2
library(Seurat)
library(ggplot2)

# Visualize the results
DimPlot(integrated_subset, reduction = "umap") + scale_x_reverse()
# extract cell names for C3
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")
cells_C3 <- WhichCells(seurat_integrated, idents = "3")

# highlight just C3
DimPlot(integrated_subset, reduction = "umap", pt.size = 1.5, cells.highlight = cells_C3,
        cols.highlight = "#0247FE", cols = "#66B032") + scale_x_reverse()+NoLegend()+NoAxes()





#txd12
FeaturePlot(integrated_subset, "HyS0042.111", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#Ncol1
FeaturePlot(integrated_subset, "HyS0008.263", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

#Laminin subunit
FeaturePlot(integrated_subset, "HyS0032.220", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#NEMCILA
FeaturePlot(integrated_subset, "HyS0030.203", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

