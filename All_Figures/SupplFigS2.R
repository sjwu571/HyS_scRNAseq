library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots
library(colorspace)



#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_v2.RData")
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



#############Dotplots##############
DefaultAssay(seurat_integrated) <- "integrated"

features <- c("HyS0010.323", "HyS0034.215", "HyS0002.407", "HyS0030.52",
              "HyS0004.64", "HyS0078.51", "HyS0020.43", "HyS0103.16",
              "HyS0030.203", "HyS0010.152", "HyS0010.9", "HyS0019.243",
              "HyS0042.80", "HyS0008.263", "HyS0002.281", "HyS0034.100",
              "HyS0056.50", "HyS0320.9", "HyS0017.39", "HyS0050.7",
              "HyS0059.86", "HyS0061.60", "HyS0034.90", "HyS0085.53",
              "HyS0013.338", "HyS0009.155", "HyS0013.62", "HyS0015.116",
              "HyS0004.446", "HyS0004.396", "HyS0041.99", "HyS0001.421",
              "HyS0014.353", "HyS0045.75", "HyS0029.183")

# Create the DotPlot
p <- DotPlot(
  seurat_integrated,
  features = features,
  scale = TRUE,
  dot.min = 0.1,
  dot.scale = 10
)

# Remove axis text but keep ticks, and add gridlines
p <- p +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica", face = "bold"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_line(color = "grey90"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)  # top, right, bottom, left
  )

p



######################UMAPs######################


DefaultAssay(seurat_integrated) <- "RNA"


##Pitx [t5275aep/XM_002164950.3] (HyS0010.323)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0010.323",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Fat4-like [CRX73260] (HyS0034.215)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0034.215",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Notch2 (Crumbs-like) [T2MDK9] (HyS0002.407)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0002.407",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Hemicentin1 like1 [CRX73261] (HyS0030.52)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0030.52",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##FiCol fibrillar collagen [Q8MUF5] (HyS0004.64)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0004.64",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Astacin3 (HyS0078.51)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0078.51",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Innexin 5 [CRX73275] (HyS0020.43)##
FeaturePlot(seurat_integrated, slot = "data", "HyS0020.43",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##frizzled 3 (HyS0103.16)
FeaturePlot(seurat_integrated, slot = "data", "HyS0103.16",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Nematocilin A (HyS0030.203)
FeaturePlot(seurat_integrated, slot = "data", "HyS0030.203",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##MORN repeat-containing protein 4-like [t12001aep] (HyS0010.152)
FeaturePlot(seurat_integrated, slot = "data", "HyS0010.152",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##adhesion G protein-coupled receptor [t9581aep] (HyS0010.9)
FeaturePlot(seurat_integrated, slot = "data", "HyS0010.9",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##calmodulin-A-like [t12424aep] (HyS0019.243)
FeaturePlot(seurat_integrated, slot = "data", "HyS0019.243",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##ARSTNd2-like (HyS0042.80)
FeaturePlot(seurat_integrated, slot = "data", "HyS0042.80",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Ncol1 (HyS0008.263)
FeaturePlot(seurat_integrated, slot = "data", "HyS0008.263",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Nematogalectin [HVAEP1.G028385] (HyS0002.281)
FeaturePlot(seurat_integrated, slot = "data", "HyS0002.281",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##NOWA (HyS0034.100)
FeaturePlot(seurat_integrated, slot = "data", "HyS0034.100",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Dkk3 (HyS0056.50, HyS0320.9, HyS0017.39)
FeaturePlot(seurat_integrated, slot = "data", "HyS0056.50",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Dkk3 (HyS0056.50, HyS0320.9, HyS0017.39)
FeaturePlot(seurat_integrated, slot = "data", "HyS0320.9",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Dkk3 (HyS0056.50, HyS0320.9, HyS0017.39)
FeaturePlot(seurat_integrated, slot = "data", "HyS0017.39",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Piwi1 (HyS0050.7)
FeaturePlot(seurat_integrated, slot = "data", "HyS0050.7",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##GNL3 (HyS0059.86)
FeaturePlot(seurat_integrated, slot = "data", "HyS0059.86",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##PCNA (HyS0061.60)
FeaturePlot(seurat_integrated, slot = "data", "HyS0061.60",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Neurocalcin (HyS0034.90)
FeaturePlot(seurat_integrated, slot = "data", "HyS0034.90",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Elav (HyS0085.53)
FeaturePlot(seurat_integrated, slot = "data", "HyS0085.53",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##RFamide (HyS0013.338)
FeaturePlot(seurat_integrated, slot = "data", "HyS0013.338",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##GLWamide (HyS0009.155)
FeaturePlot(seurat_integrated, slot = "data", "HyS0009.155",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##TSR1/hemicentin-like1 [HVAEP1.G027769] (HyS0013.62)
FeaturePlot(seurat_integrated, slot = "data", "HyS0013.62",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##mucin2-like [G010426] (HyS0015.116)
FeaturePlot(seurat_integrated, slot = "data", "HyS0015.116",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##mucin5AC (HyS0004.446)
FeaturePlot(seurat_integrated, slot = "data", "HyS0004.446",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##Rhamnospondin (HyS0004.396)
FeaturePlot(seurat_integrated, slot = "data", "HyS0004.396",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))


##Chit1 (HyS0041.99)
FeaturePlot(seurat_integrated, slot = "data", "HyS0041.99",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Dickkopf 1/2/4 A [HVAEP1.G001965] (HyS0001.421)
FeaturePlot(seurat_integrated, slot = "data", "HyS0001.421",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##matrilysin-like (t32151) (HyS0014.353)
FeaturePlot(seurat_integrated, slot = "data", "HyS0014.353",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##Irf-like (HyS0045.75)
FeaturePlot(seurat_integrated, slot = "data", "HyS0045.75",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

##allorecognition 1-like (HyS0029.183)
FeaturePlot(seurat_integrated, slot = "data", "HyS0029.183",
            order = T, min.cutoff = "q5", max.cutoff = "q99", pt.size =2, label = FALSE) + labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))



#########