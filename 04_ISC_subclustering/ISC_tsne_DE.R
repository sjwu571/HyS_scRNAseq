#' ISC subclustering DE testing
#' usually the default if used
# default Wilcox
library(Seurat)
setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/i-cell_subclustering/DE")
load("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/i-cell_subclustering/ISC_subcluter_Seurat_V6_tnse.RData")

 DefaultAssay(seurat_integrated) <- "RNA"
# DE 
de.markers <- FindAllMarkers(seurat_t, only.pos = T, 
                             min.pct = 0.5, logfc.threshold = 1, test.use = "bimod") 

# read in gene annotation files
annotation <- read.csv("hsym_combined_annotations_2.csv")
head(annotation)

# Need to reformat first column
library(stringr)
library(dplyr)

Iwant <- str_replace(annotation$gene_id, "Hsym[|]", "")
head(Iwant)
annotation$gene_id <- Iwant
head(annotation)

#de.markers <- rownames_to_column(de.markers, var = "gene_id")
head(de.markers)
head(annotation)

#tmp <- left_join(de.markers, annotation, by = "gene_id") # de.markers column 7 change name
colnames(de.markers)[7] <- "gene_id"
colnames(de.markers)

# Now it works
tmp <- left_join(de.markers, annotation, by = "gene_id") 

head(tmp)
#View(tmp)

# save
getwd()

write.csv(tmp, file = "ISC_tsne_DE_markers_min_pct_0.5_FC_1_final.csv")

# 
# # PCNA (HyS0061.60), NR2E1 (HyS0022.107), Histone H2A.Z  (HyS0044.71), MCM7B (HyS0009.219)
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0061.60","HyS0022.107", "HyS0044.71","HyS0009.219"), reduction = "tsne",
#             order = T,min.cutoff = "q5", max.cutoff = "q95")
# 
# # Piwi (HyS0050.7), nucleoplasmin-like protein ANO39 (HyS0042.86), NOP58 (HyS0155.10), HyS0073.68 (NOP56), prohibitin (HyS0017.93), GNL3 (HyS0059.86)
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0050.7","HyS0042.86", "HyS0155.10","HyS0073.68"), reduction = "tsne",
#             order = T,min.cutoff = "q5", max.cutoff = "q95")
# 
# #Sox22 (HyS0012.308), ELAV (HyS0085.53), Hippocalcin-like (HyS0034.90), NDF1 Neurogenic differentiation factor 1 (HyS0028.205)
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0012.308","HyS0085.53", "HyS0034.90","HyS0028.205"), reduction = "tsne",
#             order = T,min.cutoff = "q5", max.cutoff = "q95")
# 
# # HyS0030.178, HyS0117.11, HyS0002.164, HyS0002.228
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0030.178","HyS0117.11", "HyS0002.164","HyS0002.228"), reduction = "tsne",
#             order = T,min.cutoff = "q5", max.cutoff = "q95")
# 
# # HyS0010.86, Myc (HyS0005.84), HyS0006.58 non muscle actin
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0005.84","HyS0006.58"), reduction = "tsne",
#             order = T,min.cutoff = "q15", max.cutoff = "q95")
# 
# # C5, arrestin domain protein, Ash
# FeaturePlot(seurat_t, slot = "data", features = c("HyS0042.80", "HyS0005.437"), reduction = "tsne",
#             order = T,min.cutoff = "q5", max.cutoff = "q95")



