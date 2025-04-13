########################################################
# Cnidogenesis trajectories using Monocle3 & cytoTRACE #
########################################################

# Justin Waletich, March 2025
# justin_waletich@ufl.edu


# Loading Required Programs and Data
####################################

library(monocle3)
library(Seurat) 
library(cowplot) 
library(dplyr) 
library(ggplot2) 
library(CytoTRACE)
library(sctransform)
library(pheatmap)

############################################################################################################

# save(integrated_subset, file = "integrated_subset.RData")
load("~/Desktop/Dissertation/3_DataChapter 3/Updated_scAtlas/0_datasets/integrated_subset.RData")

# Visualize the results
DimPlot(integrated_subset, reduction = "umap", group.by = "ident")

############################################################################################################
################## Default Monocle3 pipeline with seurat object "integrated_subset"
## Uses the seurat object "integrated_subset" PCA&UMAP embeddings

# Creating object from integrated_subset to be imported into Monocle3
expression_data <- GetAssayData(integrated_subset, assay = "RNA", slot = "counts")

# Convert Seurat integrated_subset metadata to a DataFrame
cell_metadata <- integrated_subset@meta.data

# Convert Seurat integrated_subset variable features to a DataFrame
gene_annotation <- data.frame(gene_short_name = rownames(integrated_subset), row.names = rownames(integrated_subset))

# Make sure that gene_annotation and expression_data have the same number of rows (important to get Monocle3 to read dataframe)
nrow(gene_annotation)
nrow(cell_metadata)
nrow(expression_data)
common_genes <- intersect(rownames(gene_annotation), rownames(expression_data))
expression_data <- expression_data[common_genes, ]

# Create a Monocle CellDataSet object
cds <- new_cell_data_set(expression_data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Use Seurat's PCA and UMAP embeddings in Monocle
reducedDim(cds, "PCA") <- Embeddings(integrated_subset, reduction = "pca")
reducedDim(cds, "UMAP") <- Embeddings(integrated_subset, reduction = "umap")

# Cluster the cells
cds <- cluster_cells(cds)

# Learn graph function
cds <- learn_graph(cds, use_partition = F, close_loop=T)

# Monocle3 plot of data
plot_cells(cds)

# For pseudotime analysis
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


############################################################################################################
# Heatmap plot

genes_of_interest <- c("HyS0098.34", "HyS0042.111", "HyS0020.22",
                       "HyS0008.263", "HyS0020.163", "HyS0056.50",
                       "HyS0032.220", "HyS0001.466", "HyS0042.80", "HyS0030.203", "HyS0002.425",
                       "HyS0027.82")

new_row_names <- c("Pter", "Txd12", "Fkbp14-like",
                   "Ncol1", "PaxA", "Dkk3",
                   "Laminin subunit", "ZIP3-like", "Arstn2d-like", "Nematocilin A", "HyS0002.425",
                   "SCO-spondin-like")

desired_order <- c("Pter", "Txd12", "Fkbp14-like",
                   "PaxA", "Ncol1", "Dkk3",
                   "Laminin subunit", "Arstn2d-like", "ZIP3-like", "HyS0002.425", "Nematocilin A",
                   "SCO-spondin-like")

# Assumes you have already calculated pseudotime and have a Monocle3 object (cds)
pseudotime_values <- pseudotime(cds)

# Get the gene expression matrix for your genes of interest
expr_matrix <- as.matrix(exprs(cds)[genes_of_interest, ]) # extracts expression matrix data for genes of interest
sorted_cells <- names(sort(pseudotime_values)) # creates a vector that orders cells from smallest to largest pseudotime value
expr_matrix_sorted <- expr_matrix[, sorted_cells] # sorts our expression matrix and orders it from lowest pseudotime cells to largest
expr_matrix_scaled <- t(scale(t(expr_matrix_sorted))) # gets the z-score expression of each gene in every cell (normalizes the data)

# Set max and min limits for scaling
max_value <- 1
min_value <- 0

# Assign the new row names to your matrix
rownames(expr_matrix_scaled) <- new_row_names

# Reorder the rows of the matrix according to your custom order
expr_matrix_ordered <- expr_matrix_scaled[desired_order, ]

# Define the color breaks based on the range you want
breaks <- seq(min_value, max_value, length.out = 16)

# Create the heatmap
pheatmap(
  expr_matrix_ordered, 
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE, 
  show_rownames = FALSE,
  color = colorRampPalette(c("#3B014A", "#b73779", "#f46d43", "#fee08b"))(length(breaks) - 1),  # Color scale
  breaks = breaks,  # Set the breaks for color scaling
  cellheight = 12.5
)

############################################################################################################
### cytoTRACE

expression_matrix1 <- as.matrix(integrated_subset@assays$RNA@counts) # Get counts matrix from Seurat object

expression_df1 <- as.data.frame(expression_matrix1) # Make it a dataframe

results1 <- CytoTRACE(expression_df1) # Run CytoTRACE

plotCytoTRACE(results1) # Plot the results using CytoTRACE's clustering

# Extract the CytoTRACE scores
cytotrace_scores1 <- results1$CytoTRACE

# Combine the CytoTRACE scores
combined_scores <- c(cytotrace_scores1) 
combined_scores <- combined_scores[colnames(integrated_subset)] # Add CytoTRACE scores for each cell to Seurat object

# Add CytoTRACE scores to the original Seurat object as metadata
integrated_subset$CytoTRACE <- combined_scores

# Reverse CytoTRACE scores so that high means more differentiated and low means less differentiated
integrated_subset[['differentiation_scores']] <- 1 - integrated_subset[['CytoTRACE']] 


# Code for CytoTRACE plot 
FeaturePlot(integrated_subset, features = "differentiation_scores", reduction = "umap", pt.size = 1.5, order = TRUE) + 
  theme(plot.title = element_text(face = "bold.italic", size = 35)) +
  labs(title = "cytoTRACE") +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),    # Remove Y axis text
    axis.line = element_blank(),      # Remove axis lines
    axis.ticks = element_blank()      # Remove axis ticks
  ) + 
  scale_color_gradientn(colors = rev(c("#902843", "#B4404D", "#DC7449", "#EAAB64", "#F6DE87", "#F4D835", "#E7F198", "#B9DCA5", "#82BEA6", "#536FAF", "#5A57A3"))) +
  theme(
    legend.title = element_text(size = 12, face = "bold"),  # Customize legend title
    legend.key = element_blank(),  # Remove the key background
    legend.key.width = unit(1, "cm"),  # Adjust the key width size if needed
    legend.key.height = unit(0.4, "cm"),  # Adjust the key height size if needed
    legend.text = element_blank(),  # Remove legend text (numbers)
    legend.ticks = element_blank(),  # Remove ticks
    legend.position = c(0.6, 0.03),  # Position legend at the bottom right
    legend.direction = "horizontal",  # Make the legend horizontal
    legend.box = "horizontal"  # Arrange legend items in a horizontal box
  )

