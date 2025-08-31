###Figure 5 v4_3 (use the RNA Assay for Featureplots, maxcutoff q99)
###Use hex codes for colors#####
###export image as eps file, 1212 x 1000####

library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots

#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final.RData")
# see UMAP
DimPlot(seurat_integrated, label = F, label.size = 5, pt.size=1)
# need to run this to change the default assay for the current UMAP
DefaultAssay(seurat_integrated) <- "RNA"

tmp <- seurat_integrated
table(Idents(tmp))
levels(x = tmp)

cluster.ids <- c("epithelial cells", "nematocytes", "nematocytes", "i-cells", "nematocytes",
                 "epithelial cells", "nematoblasts", "neurons", "nematocytes", "immune cells",
                 "gland cells", "gland cells", "gland cells", "gland cells", "neurons",
                 "nematoblasts", "nematoblasts", "epithelial cells", "epithelial cells") 
names(cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, cluster.ids)

a <- c('#e03672', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey','lightgrey', 'lightgrey')

DimPlot(tmp, reduction = "umap", order = TRUE, pt.size = 2, cols = a) + NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#Astacin3_HyS0078.51_magenta
FeaturePlot(seurat_integrated, "HyS0078.51", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#Fat1_HyS0048.57_yellow
FeaturePlot(seurat_integrated, "HyS0048.57", order = TRUE, cols = c("lightgrey", "yellow"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 


#HyS0001.363_magenta
FeaturePlot(seurat_integrated, "HyS0001.363", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 
