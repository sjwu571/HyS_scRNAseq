### Cite: https://academic.oup.com/bioinformatics/article/37/19/3374/6198103 ####
### https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html ###
# deploy seurat object on shiny.io ###

# reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
#            "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
# newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
# if(length(newPkg)){install.packages(newPkg)}
# 
# # If you are using h5ad file as input, run the code below as well
# # reticulate::py_install("anndata")
# 
 devtools::install_github("SGDDNB/ShinyCell", force = T)

library(Seurat)
library(ShinyCell)
library(rsconnect)

#getExampleData()                      
# seurat_integrated = readRDS("readySeu_rset.rds")
load("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/Seurat_integrated_features3k_PC50_dims20_res0.5_final_clean.RData")
DefaultAssay(seurat_integrated) <- "RNA" 
seurat_integrated <- FindVariableFeatures(seurat_integrated)

# integrated assay only has 3000 features/genes

DimPlot(seurat_integrated, label = T)
getwd()
# set working directory
setwd("~/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/ShinyCell/V5")

scConf = createConfig(seurat_integrated,
                      meta.to.include = c("integrated_snn_res.0.3", "sample")) 

# new folder "ShinyApp" created in the current directory
makeShinyApp(seurat_integrated, scConf, gex.assay = "RNA", gene.mapping = FALSE,
             shiny.title = "Updated Hsym single-cell atlas",
             default.gene1="HyS0008.263", default.gene2="HyS0030.203",
             default.multigene=c("HyS0042.111", "HyS0008.263", "HyS0032.220", "HyS0030.203")
             ) 


# this step fails if assay  = SCT. gex.assay has to be either "RNA" or "integrated"
# Warning message:
#In makeShinyFiles(obj = obj, scConf = scConf, gex.assay = gex.assay[1],  :
 #                   Variable genes for seurat object not found! Have you ran `FindVariableFeatures` or `SCTransform`?

#install.packages('rsconnect')
#You need a ShinyApp account for the token and secrets
rsconnect::setAccountInfo(name='sjwu571', token='', secret='')

rsconnect::deployApp('/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/ShinyCell/V5/shinyApp',
                     appName = "hys-umap")


# https://sjwu571.shinyapps.io/hys-umap/
