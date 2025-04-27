#' ISC subclustering DE testing
#' usually the default if used
# default Wilcox
library(Seurat)
setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/i-cell_subclustering/DE")
load("../ISC_subcluter_Seurat_V6_tnse.RData")

DefaultAssay(seurat_t) <- "RNA"

# DE 
de.markers <- FindAllMarkers(seurat_t, only.pos = T, 
                             min.pct = 0.3, logfc.threshold = 1) 

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

write.csv(tmp, file = "ISC_tsne_DE_markers_min_pct_0.3_FC_1_final.csv")

