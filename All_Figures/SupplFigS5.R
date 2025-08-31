# Load your Seurat object
load("~/Desktop/R files/sub_cluster_7_14_resolution0.1")

#load Jing-wei's UMAP saved on desktop
load("~/Desktop/R files/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")

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

a <- c('lightgrey', 'lightgrey', 'lightgrey', 'lightgrey', '#9D2EFF','lightgrey', 'lightgrey')

DimPlot(tmp, reduction = "umap", order = TRUE, pt.size = 2, cols = a) + NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 

##################nerual subcluster Dimplot####################################

DimPlot(sub_cluster7_14_resolution0.1, order = TRUE, pt.size = 2) + NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))

################### Visualize ORIGINAL clusters (7 vs 14)#####################

DimPlot(sub_cluster7_14_resolution0.1, group.by = "orig_cluster", pt.size = 2) + NoLegend() +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) + ggtitle(NULL)

#########################Neuropeptide gene expression##########################

DefaultAssay(sub_cluster7_14_resolution0.2) <- "RNA"

###

#HyS0009.155_GLWamide
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0009.155", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0013.338_RFamide
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0013.338", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0018.64
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0018.64", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0019.179
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0019.179", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0052.141
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0052.141", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0078.10
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0078.10", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0044.153
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0044.153", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0049.55
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0049.55", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0055.9
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0055.9", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))


#HyS0051.139
FeaturePlot(sub_cluster7_14_resolution0.2, "HyS0051.139", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0028.162
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0028.162", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

#HyS0062.54
FeaturePlot(sub_cluster7_14_resolution0.1, "HyS0062.54", order = TRUE, 
            min.cutoff = "q5", max.cutoff = "q99", slot = "data", pt.size = 2, label = TRUE) + 
  theme(axis.text.x = element_text(face = "bold", size = 20, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 20, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.25, "cm"))

