###Figure 1 v7_3 (changing colors of clusters to make them brighter/bolder, changed neuron cluster color, also remove axes from main atlas UMAP, use the RNA Assay for Featureplots, maxcutoff q95)
###Use hex codes for colors#####
###export image as eps file, 1212 x 1000####

library(Seurat)
library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots
library(colorspace)


#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")
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

a <- c('#e03672', '#66B032', '#0247FE', '#66B032', '#8601AF','#EBEB00', '#F28522')

DimPlot(tmp, reduction = "umap", label.size = 5, pt.size = 1.5, cols = a) + 
  NoLegend() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "cm")
  )


#Ncol1_HyS0008.263_green
original_color <- "#66B032"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0008.263", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    

#NematocilinA_HyS0030.203_Green
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#66B032"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0030.203", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    


#Elav_HyS0085.53_purple
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#9D2EFF"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0085.53", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

#Mucin-5AC_HyS0004.446_orange
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#F28522"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0004.446", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    


#Chit1_HyS0041.99_orange
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#F28522"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0041.99", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    


#Fat1_HyS0048.57_Red
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#e03672"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0048.57", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    



#Astacin3_HyS0078.51_Red
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#e03672"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0078.51", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    



#MCM5_HyS0033.31_Blue
DefaultAssay(seurat_integrated) <- "RNA"
original_color <- "#0247FE"
darker_color <- darken(original_color, amount = 0.2)

FeaturePlot(seurat_integrated, "HyS0033.31", order = TRUE, cols = c("lightgrey", darker_color), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))    


