###export image as eps file, 1212 x 1000####
###Use RNA assay as default, maxcutoff q99#####

##Plot genes that show gradient expression in C0###
library(Seurat) # main program needed
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

##making everything grey except c5###
cluster.ids <- c("epithelial cells0", "nematocytes1", "nematocytes2", "i-cells3", "nematocytes4",
                 "epithelial cells5", "nematoblasts6", "neurons7", "nematocytes8", "immune cells9",
                 "gland cells10", "gland cells11", "gland cells12", "gland cells13", "neurons14",
                 "nematoblasts15", "nematoblasts16", "epithelial cells17", "epithelial cells18") 
names(cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, cluster.ids)

a <- c('lightgrey','lightgrey','lightgrey','lightgrey','lightgrey','#e03672', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey','lightgrey', 'lightgrey','lightgrey','lightgrey','lightgrey','lightgrey','lightgrey','lightgrey','lightgrey')

DimPlot(tmp, reduction = "umap", label.size = 5, pt.size = 1.0, cols = a) + 
  NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0006.325
FeaturePlot(seurat_integrated, slot = "data", "HyS0006.325",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0028.60
FeaturePlot(seurat_integrated, slot = "data", "HyS0028.60",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))





##########################################################################
##making everything grey except c0###
cluster.ids <- c("epithelial cells", "nematocytes", "nematocytes", "i-cells", "nematocytes",
                 "nematocytes", "nematoblasts", "neurons", "nematocytes", "immune cells",
                 "gland cells", "gland cells", "gland cells", "gland cells", "neurons",
                 "nematoblasts", "nematoblasts", "nematocytes", "nematocytes") 
names(cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, cluster.ids)

a <- c('#e03672', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey','lightgrey', 'lightgrey')

DimPlot(tmp, reduction = "umap", label.size = 5, pt.size = 1.0, cols = a) + 
  NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))



#HyS0008.375
FeaturePlot(seurat_integrated, slot = "data", "HyS0008.375",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0021.141
FeaturePlot(seurat_integrated, slot = "data", "HyS0021.141",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0054.100
FeaturePlot(seurat_integrated, slot = "data", "HyS0054.100",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0034.215
FeaturePlot(seurat_integrated, slot = "data", "HyS0034.215",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0053.75
FeaturePlot(seurat_integrated, slot = "data", "HyS0053.75",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0087.37
FeaturePlot(seurat_integrated, slot = "data", "HyS0087.37",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0244.2
FeaturePlot(seurat_integrated, slot = "data", "HyS0244.2",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))
