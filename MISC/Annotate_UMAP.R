# Jingwei, June 8, 2023
# rename clusters in UMAP
# updated 12/2/23
# 1/21/24 Fig S1, detailed annotation of each cluster
# 9/29/24 cluster names update


library(Seurat)
load("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_v2.RData")

#' from Seurat 
#' new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#'                     "NK", "DC", "Platelet")
#'names(new.cluster.ids) <- levels(pbmc)
#'pbmc <- RenameIdents(pbmc, new.cluster.ids)
#'DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#+ code block
#+ https://docs.google.com/spreadsheets/d/1rAb6J75VfmRHpRLpm195eNQKcgrEfVV1UaPwaRnI9YU/edit#gid=1284240333
# From C0 to C18
new.cluster.ids <- c("Endodermal epithelial", "Mature desmoneme", "Developing desmoneme", "I-cells/progenitors", #C0-C3
                     "Developing nematocytes", "Ectodermal epithelial", "Nematoblasts", "Ganglion neurons", #C4-C7
                     "Mature eurytele", "Unknown", "Mucous GCs", "Mucous GCs", #C8 - 11
                      "Zymogen GCs", "Zymogen GCs", "Sensory neurons", "Nematoblasts", #C12-C15
                     "Nematoblasts","Endodermal epithelial","Stolon-specific Epi") #C16-C18
# Remove C19
Iwant <- subset(seurat_integrated, idents = 19, invert = T)
DimPlot(Iwant)
seurat_integrated <- Iwant
                     
names(new.cluster.ids) <- levels(seurat_integrated)
#View(seurat_integrated@meta.data)
dim(seurat_integrated@meta.data)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
DimPlot(seurat_integrated, label = T, pt.size = 0.5, label.size = 5) + NoLegend()

#+ 
#+ Change colors
#+ https://stackoverflow.com/questions/63867603/how-to-change-the-default-color-scheme-of-seurat-dimplot
my_cols <- c('Endodermal epithelial'='#FF1F5B', 'Ectodermal epithelial'='#FF1F5B',
             'Mature desmoneme'='#00CD6C', 'Developing desmoneme'='#00CD6C', 'Developing nematocytes'='#00CD6C', 'Nematoblasts'='#00CD6C', 'Mature eurytele' = '#00CD6C',
               'I-cells/progenitors'='#009ADE',
               'Ganglion neurons'='#AF58BA','Sensory neurons' ='#AF58BA',
               'Mucous GCs'='#FFC61E', 'Zymogen GCs'='#FFC61E',
             'Unknown'='#F28522', 'Stolon-specific Epi'='#A6761D')

# DimPlot(seurat_integrated,
#         cols = my_cols, label=TRUE , repel=TRUE, pt.size = 0.5, label.size = 5) + NoLegend()

DimPlot(seurat_integrated,
        cols = my_cols, label=F) + NoLegend()





# ###### 1/21/24 ########
# 
# library(Seurat)
# load("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_v2.RData")
# 
# DimPlot(seurat_integrated)
# 
# # Remove C19
# Iwant <- subset(seurat_integrated, idents = 19, invert = T)
# 
# DimPlot(Iwant, label = T, pt.size = 0.5, label.size = 5) + NoAxes()
# 
# # new.cluster.ids <- c("Endodermal epithelial", "Mature desmoneme", "Developing desmoneme", "I-cells/progenitors", #C0-C3
# #                      "Developing nematocytes", "Epithelial", "Nematoblasts", "Neuron Type A", #C4-C7
# #                      "Mature eurytele", "Unknown", "Mucous GCs", "Mucous GCs", #C8 - 11
# #                      "Zymogen GCs", "Zymogen GCs", "Neuron Type B", "Nematoblasts", #C12-C15
# #                      "Nematoblasts","Endodermal epithelial","Unknown") #C16-C18
# # 



# Supp fig. 
# revert annotations to numbers
load("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_v2.RData")
Iwant <- subset(seurat_integrated, idents = 19, invert = T)
DimPlot(Iwant)
seurat_integrated <- Iwant

my_cols <- c('Endodermal epithelial'='#FF1F5B', 'Ectodermal epithelial'='#FF1F5B',
             'Mature desmoneme'='#00CD6C', 'Developing desmoneme'='#00CD6C', 'Developing nematocytes'='#00CD6C', 'Nematoblasts'='#00CD6C', 'Mature eurytele' = '#00CD6C',
             'I-cells/progenitors'='#009ADE',
             'Ganglion neurons'='#AF58BA','Sensory neurons' ='#AF58BA',
             'Mucous GCs'='#FFC61E', 'Zymogen GCs'='#FFC61E',
             'Unknown'='#F28522', 'Stolon-specific Epi'='#A6761D')


  my_cols <- c(
    '0' = '#FF1F5B',
    '1' = '#00CD6C',
    '2' = '#00CD6C',
    '3' = '#009ADE',
    '4' = '#00CD6C',
    '5' = '#FF1F5B',
    '6' = '#00CD6C',
    '7' = '#AF58BA',
    '8' = '#3498DB',
    '9' = '#F28522',
    '10' = '#FFC61E',
    '11' = '#FFC61E',
    '12' = '#FFC61E',
    '13' = '#FFC61E',
    '14' = '#AF58BA',
    '15' = '#00CD6C',
    '16' = '#00CD6C',
    '17' = '#FF1F5B',
    '18' = '#A6761D'
  )


# DimPlot(seurat_integrated,
#         cols = my_cols, label=TRUE , repel=TRUE, pt.size = 0.5, label.size = 5) + NoLegend()

DimPlot(seurat_integrated,
        cols = my_cols, label=T, label.size = 8) + NoLegend()



