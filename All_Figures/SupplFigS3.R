#load packages
library(Seurat) # main program needed
library(cowplot) # used if making multiplots
library(dplyr) # used in various functions
library(ggplot2) # used to customize plots


#load Cnidogenesis UMAP saved on desktop
load("~/Desktop/R files/integrated_subset.RData")
# need to run this to change the default assay for the current UMAP
DefaultAssay(integrated_subset) <- "RNA"

#Fkbp14
FeaturePlot(integrated_subset, "HyS0020.22", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  
#Dkk3
FeaturePlot(integrated_subset, "HyS0056.50", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

#Arstnd2
FeaturePlot(integrated_subset, "HyS0042.80", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm")) 
#c1 marker
FeaturePlot(integrated_subset, "HyS0002.425", order = TRUE, cols = c("lightgrey", "magenta"), 
            min.cutoff = "q5",max.cutoff ="q99", slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  

#c8 marker_yellow
FeaturePlot(integrated_subset, "HyS0027.82", order = TRUE, cols = c("lightgrey", "yellow"), 
            min.cutoff = "q5", max.cutoff ="q99",slot = "data", pt.size = 2, label = FALSE) + 
  scale_x_reverse() +  # <-- this flips the plot horizontally
  labs(title = NULL) +
  theme(axis.text.x = element_text(face = "bold", size = 30, vjust = 0),
        axis.text.y = element_text(face = "bold", size = 30, hjust = 1),
        axis.title.y = element_text(face = "bold", size = 30),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 2.5),
        axis.ticks = element_line(linewidth = 2.5),
        axis.ticks.length = unit(0.25, "cm"))  
