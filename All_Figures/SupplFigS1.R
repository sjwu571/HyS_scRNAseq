###########UMAP##############
library(Seurat)
library(ggplot2)

DimPlot(seurat_integrated, label = FALSE, label.size = 5, pt.size = 1) +
  theme(
    axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
    axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
    axis.title.y = element_text(face = "bold", size = 30),
    axis.title.x = element_text(face = "bold", size = 30),
    axis.line = element_line(linewidth = 2.5),
    axis.ticks = element_line(linewidth = 2.5),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = "none"
  )


####Give me Hex color codes for the clusters####
# Make the plot
p <- DimPlot(seurat_integrated, label = FALSE, pt.size = 1)

# Extract the ggplot build object
pb <- ggplot_build(p)

# This contains the color scale mapping
color_scale <- pb$plot$scales$get_scales("colour")

# Levels of the clusters
cluster_ids <- levels(Idents(seurat_integrated))

# Hex codes in correct order
hex_codes <- color_scale$palette(length(cluster_ids))

# Combine into a data frame
color_map <- data.frame(cluster = cluster_ids, color = hex_codes)
color_map


##Dotplot_no axes titles
library(Seurat)
library(ggplot2)

DefaultAssay(seurat_integrated) <- "integrated"

features <- c("HyS0008.263","HyS0030.203","HyS0085.53",
              "HyS0004.446","HyS0041.99","HyS0048.57",
              "HyS0078.51","HyS0033.31")

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
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_line(color = "grey90"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)  # top, right, bottom, left
  )

p







