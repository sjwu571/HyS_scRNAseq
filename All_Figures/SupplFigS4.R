###Supplementary FigS4 v2_1 (use the RNA Assay for Featureplots, maxcutoff q99, use integrated assay for dotplots)
###Use hex codes for colors#####



###export image as eps file, 1212 x 1000####
###Use RNA assay as default, maxcutoff q99#####
##Export as 1212 x 1000###

library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots

#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")
# see UMAP
DimPlot(seurat_integrated, label = F, label.size = 5, pt.size=1)
# need to run this to change the default assay for the current UMAP
DefaultAssay(seurat_integrated) <- "RNA"

#HyS0062.54_magenta
FeaturePlot(seurat_integrated, "HyS0062.54", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

#############Overall atlas UMAP#############
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



############dotplot############

###Supplementary FigS4 v2 (use the RNA Assay for Featureplots, maxcutoff q99, use integrated assay for dotplots)
###Use hex codes for colors#####

##dotplots

library(Seurat)
library(ggplot2)
library(rlang)  # for quo()

DefaultAssay(seurat_integrated) <- "integrated"

features <- c("HyS0009.155","HyS0013.338","HyS0018.64",
              "HyS0019.179","HyS0052.141","HyS0078.10",
              "HyS0044.153","HyS0049.55", "HyS0055.9", "HyS0051.139",
              "HyS0028.162", "HyS0062.54")

gene_colors <- c(
  "HyS0009.155" = "#be1e2d",
  "HyS0013.338" = "#be1e2d",
  "HyS0018.64"  = "#2e3192",
  "HyS0019.179" = "#be1e2d",
  "HyS0052.141" = "#be1e2d",
  "HyS0078.10"  = "#be1e2d",
  "HyS0044.153" = "#9D2EFF",
  "HyS0049.55"  = "#2e3192",
  "HyS0055.9"   = "#be1e2d",
  "HyS0051.139" = "#9D2EFF",
  "HyS0028.162" = "#9D2EFF",
  "HyS0062.54"  = "#2e3192"
)

# Build the default DotPlot
p <- DotPlot(
  seurat_integrated,
  features = features,
  scale = TRUE, dot.min = 0.1, dot.scale = 15
)

# Ensure x-axis follows your input order
p$data$features.plot <- factor(p$data$features.plot, levels = features)

# Re-map color from continuous avg.exp.scaled -> discrete features.plot
p$layers[[1]]$mapping$colour <- quo(features.plot)

# Apply custom colors and theme
p <- p +
  scale_color_manual(values = gene_colors, name = "Gene") +
  guides(size = guide_legend(title = "% Expressing")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica", face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, family = "Helvetica", face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_line(color = "grey90"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  ) +
  labs(y = "cluster", x = NULL)  # <- make sure this is part of the chain

p


#####################UMAPS########################################
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

#GLWamide
FeaturePlot(seurat_integrated, "HyS0009.155", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#RFamide
FeaturePlot(seurat_integrated, "HyS0013.338", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#Involucrin+neurons_from Diaz
FeaturePlot(seurat_integrated, "HyS0018.64", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#uncharacterized
FeaturePlot(seurat_integrated, "HyS0019.179", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#PRXamide?
FeaturePlot(seurat_integrated, "HyS0052.141", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#Noblasthit
FeaturePlot(seurat_integrated, "HyS0078.10", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#Uncharacterized
FeaturePlot(seurat_integrated, "HyS0044.153", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#PQRFVamide/myb-like
FeaturePlot(seurat_integrated, "HyS0049.55", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#uncharacterized
FeaturePlot(seurat_integrated, "HyS0055.9", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#No blast hit
FeaturePlot(seurat_integrated, "HyS0051.139", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#Hypothetical protein
FeaturePlot(seurat_integrated, "HyS0028.162", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#No blast hit
FeaturePlot(seurat_integrated, "HyS0062.54", order = TRUE, cols = c("grey", "#8601AF"), 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))



