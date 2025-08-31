###Supplementary Fig 4 v1_3###
###export image as eps file, 1212 x 1000####
###Use RNA assay as default, maxcutoff q99#####

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

#HyS0045.75_magenta
FeaturePlot(seurat_integrated, "HyS0045.75", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#HyS0053.57_magenta
FeaturePlot(seurat_integrated, "HyS0053.57", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 



