# https://satijalab.org/seurat/reference/doheatmap

# DoHeatmap(
#   object,
#   features = NULL,
#   cells = NULL,
#   group.by = "ident",
#   group.bar = TRUE,
#   group.colors = NULL,
#   disp.min = -2.5,
#   disp.max = NULL,
#   slot = "scale.data",
#   assay = NULL,
#   label = TRUE,
#   size = 5.5,
#   hjust = 0,
#   angle = 45,
#   raster = TRUE,
#   draw.lines = TRUE,
#   lines.width = NULL,
#   group.bar.height = 0.02,
#   combine = TRUE
# )


# load seuart V2
load("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_v2.RData")
# make sure right ident
library(Seurat)
DimPlot(seurat_integrated)
# add cell types
# new.cluster.ids <- c("Endo Epi", "Desmoneme NCs", "Desmoneme NCs", "Interstitial SCs",
#                      "Desmoneme NCs", "Epi", "Desmoneme NBs", "Neural", #C4-C7
#                      "Eurytele NCs", "unknown", "Spumous Muc GCs", "Granular Muc GCs", #C8 - 11
#                      "Zymogen GCs", "Zymogen GCs", "Neural", "Eurytele NBs",
#                      "Desmoneme NBs", "Endo Epi", "Ecto Epi","19")
# names(new.cluster.ids) <- levels(seurat_integrated)
# View(seurat_integrated@meta.data)
# dim(seurat_integrated@meta.data)

# seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
# 
# DimPlot(seurat_integrated, label = T, label.size = 3) + NoLegend()


# DoHeatmap(seurat_integrated,
#           features = c("HyS0008.263","HyS0042.80","HyS0005.468"), # add more marker genes
#           group.by = "ident", # delete C19
#           size = 3, draw.lines = T) +
#   custom_colors
# 
# 
# custom_colors <- scale_fill_gradient(
#   low = "white",   # Define the color for low values
#   high = "red",   # Define the color for high values
#   breaks = c(-2, 0, 2),  # Define the breaks for the color scale
#   labels = c("Low", "Medium", "High")  # Define labels for the breaks
# )
# 

# get all clusters from nematoblasts to mature cnidocytes
Iwant <- subset(seurat_integrated, idents = c("15","16","6","4","2","1","8"))
Iwant
levels(Iwant) 
#[1] "1"  "2"  "4"  "6"  "8"  "15" "16"
DimPlot(Iwant, label = T)

# 1-7 early to late nematogenesis
new.cluster.ids <- c("6","5","4","3","7","1","2")

names(new.cluster.ids) <- levels(Iwant)
View(seurat_integrated@meta.data)
dim(seurat_integrated@meta.data)

# key in renaming levels. same function used in renaming cluster names
Iwant <- RenameIdents(Iwant, new.cluster.ids)
Idents(Iwant)

# heatmap
# library(ggplot2)
# custom_colors <- scale_fill_gradient(
#   low = "#3300CC",   # Define the color for low values
#   high = "#FF3300",   # Define the color for high values
#  # breaks = c(-2, 0, 2),  # Define the breaks for the color scale
#   breaks = c(0, 3, 6),  # Define the breaks for the color scale
#   
#   labels = c("Low", "Medium", "High")  # Define labels for the breaks
# )

# picks genes that express at different part of the trajectory from early to late
# NOWA (HyS0034.100), Dickkopf related-protein 3-like (HyS0056.50)
# Nematocyst_outer_wall_antigen (HyS0005.226), Nematoglactin B (HyS0002.32)
# Ncol1 (HyS0008.263), Nematoglactin (HyS0002.281), NR2E1 (HyS0022.107)
# laminin subunit alpha-1-like (HyS0032.220), PFB0145c-like (HyS0044.107), PFB0765w-like (HyS0111.9)
# ARSTNd2-like (HyS0042.80), Voltage-gated_potassium_channel_alpha_subunit (HyS0026.207), laminin subunit alpha-1-like (HyS32.220)
# Nematocilin A (HyS0030.203), polycystin-1-like (HyS0040.103), transmembrane_protein_180-like (HyS0014.58)
# Nematocilin A (HyS0030.203), calcium-dependent secretion activator 1-like (HyS0076.19), Potassium voltage-gated channel protein (HyS0003.504), metabotropic glutamate receptor 3-like (HyS0006.178), uncharacterized protein LOC105846444 (HyS0076.32)


##### 5/27/24 #####

levels(Iwant) <- as.factor(c(1:7))

DoHeatmap(Iwant,
          features = c("HyS0022.107","HyS0056.50","HyS0008.263",
                       "HyS0002.281","HyS0005.226","HyS0002.32",
                       "HyS0032.220","HyS0042.80",
                       "HyS0030.203","HyS0040.103", "HyS0012.244","HyS0007.266",
                       "HyS0003.504", "HyS0027.82",
                       "HyS0002.248"), 
          group.by = "ident",
          size = 3, draw.lines = T, assay = "RNA",slot = "counts") + custom_colors

# 10/6/24, modify the figure with viridis color
#group.bar = TRUE,    label = TRUE,
library(viridis)


DoHeatmap(Iwant, 
          features = c("HyS0022.107", "HyS0056.50", "HyS0008.263",
                       "HyS0002.281", "HyS0005.226", "HyS0002.32",
                       "HyS0032.220", "HyS0042.80",
                       "HyS0030.203", "HyS0040.103", "HyS0012.244", "HyS0007.266",
                       "HyS0003.504", "HyS0027.82",
                       "HyS0002.248"), 
          group.by = "ident", 
          size = 3, 
          draw.lines = TRUE, group.bar = F, label = F,
          assay = "SCT", 
          slot = "data") +
  scale_fill_viridis(option = "C",  # Choose an appropriate viridis palette (A, B, C, or D)
                     name = "Expression")  # Customize legend title


# 10 x 4 inches look good. need to remove some genes that don't look good.

DoHeatmap(Iwant, 
          features = c("HyS0056.50", "HyS0008.263",
                       "HyS0002.281", "HyS0005.226", "HyS0002.32",
                       "HyS0032.220", "HyS0042.80",
                       "HyS0030.203", "HyS0040.103", "HyS0012.244", "HyS0007.266",
                       "HyS0003.504","HyS0002.248"), 
          group.by = "ident", 
          size = 3, 
          draw.lines = TRUE, group.bar = F, label = F,
          assay = "SCT", 
          slot = "data") +
  scale_fill_viridis(option = "C",  # Choose an appropriate viridis palette (A, B, C, or D)
                     name = "Expression")  # Customize legend title


