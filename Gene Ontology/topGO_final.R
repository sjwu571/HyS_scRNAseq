### Gene Ontology Enrichment Analysis ###
# TopGO
# Follow DE step

library(dplyr)

setwd("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/topGO")

#Import PANNZER gene annotations
descriptions_original <- read_tsv("Hsym_v1.0_descriptions.out", col_names = TRUE, skip_empty_rows = TRUE)

descriptions <- descriptions_original %>% dplyr::select(qpid, desc, genename)
colnames(descriptions) <- c("gene_id", "pannzer_description", "pannzer_gene_name")
head(descriptions)
go_terms <- read_tsv("Hsym_v1.0_GO_terms.out", col_names = TRUE, skip_empty_rows = TRUE)  

head(go_terms)

go_MF <- go_terms[go_terms$ARGOT_rank == 1,] %>% filter(ontology == "MF") %>% dplyr::select(qpid, goid)
colnames(go_MF) <-  c("gene_id", "molecular_function")
go_BP <-go_terms[go_terms$ARGOT_rank == 1,] %>% filter(ontology == "BP") %>% dplyr::select(qpid, goid)
colnames(go_BP) <-  c("gene_id", "biological_processes")
go_CC <- go_terms[go_terms$ARGOT_rank == 1,] %>% filter(ontology == "CC") %>% dplyr::select(qpid, goid)
colnames(go_CC) <-  c("gene_id", "cellular_component")
annotations <- left_join(descriptions, go_MF, by = "gene_id") %>% left_join(go_BP, by = "gene_id") %>% left_join(go_CC, by = "gene_id")

head(annotations)
 #write.table(annotations, "annotations.txt", sep = '\t', quote = F)

 #BiocManager::install("topGO")
library('topGO')
library(ggplot2)

create_GO_data <- function(DE_result){
  # DE_result is a vector of gene_ids of all upregulated genes
  results <- data.frame(gene_id = DE_result$gene_id)
  # Get rid of genes with no GO term annotation
  results_with_GO <- results[results$gene_id %in% go_BP$gene_id,]
  # Make named factor that indicates which genes are interesting (upregulated) and which are not from the gene universe
  # The gene universe is all genes with GO terms from the pannzer annotation file
  
  # First, format GO terms
  go_BP_formatted <- go_BP
  go_BP_formatted$GO_BP <- paste("GO:", go_BP_formatted$biological_processes, sep = "")
  # The list of gene_ids in go_BP formatted is the gene universe
  gene_universe_names <- go_BP_formatted$gene_id
  
  # Make factor and name each element with gene name
  gene_factor <- factor(as.integer(gene_universe_names %in% results_with_GO))
  names(gene_factor) <- gene_universe_names
  
  # Write function that makes peoperly formatted table of GO term annotations  
  make_gene2GO <- function(annotation){
    gene_id <- annotation$gene_id
    GO <- annotation$GO_BP
    table <- list(c())
    for (i in 1:length(gene_id)){
      table[[i]] <- GO[i]
    }
    names(table) <- gene_id
    table
  }
  go_BP_list <-  make_gene2GO(go_BP_formatted)
  GOdata <- new("topGOdata",
                description = "TEST",
                ontology = "BP",
                allGenes = gene_factor,
                annot = annFUN.gene2GO,
                gene2GO = go_BP_list)
}



run_all_tests <- function(topGO_data_object, title = ""){
  resultFisher <- runTest(topGO_data_object, algorithm = "classic", statistic = "fisher")
  result_table <- GenTable(topGO_data_object, classicFisher = resultFisher, topNodes = 100)
  printGraph(topGO_data_object, resultFisher, firstSigNodes = 30, fn.prefix = title, useInfo = "all", pdfSW = TRUE)
  
  goEnrichment <- result_table
  goEnrichment <- goEnrichment[1:30,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  goEnrichment$classicFisher <- as.numeric(gsub("<\\s*", "", goEnrichment$classicFisher))
  #goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
  goEnrichment$classicFisher
  #file(result_table) # unused argument path = paste(title, "_table.tsv", sep = "")
  
  # write_tsv is deprecated. use file instead
  
  plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab("Biological process") +
    ylab("Enrichment") +
    ggtitle(title) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$classicFisher)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=10, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=9, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=9, face="bold", vjust=0.5),
      axis.title=element_text(size=10, face="bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=9),  #Text size
      title=element_text(size=9)) +
    guides(colour=guide_legend(override.aes=list(size=1))) +
    coord_flip()
  
  pdfname <- paste(title, "_graph.pdf", sep ="" )
  pdf(pdfname)
  print(plot)
  dev.off()
  
  list(result_table = result_table, plot = plot)
  write.table(result_table, "result_table", sep = "\t", quote = F)
}


### one cluster at a time ###
# DE_result <- read.table("/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/DE/C_0_DE_gene_ID.txt", header = T)
# GOdata <- create_GO_data(DE_result)
# run_all_tests(GOdata, "C0")
# 

#### loop through all 19 clusters ####
# DEGs gene_ID from the large DE table

# Define the base directory path
base_dir <- "/Users/jsong/OneDrive - University of Florida/Schnitzler_Lab/Jingwei/Whitney/Manuscript/10X_/CCA integration/DE/"

# Loop through clusters 0 to 18
for (i in 0:18) {
  # Construct the file name
  file_name <- paste0("C_", i, "_DE_gene_ID.txt")
  full_path <- paste0(base_dir, file_name)
  
  # Read the data
  DE_result <- read.table(full_path, header = TRUE)
  
  # Create GO data
  GOdata <- create_GO_data(DE_result)
  
  # Run tests with the appropriate cluster name
  # Using gsub to remove underscore for the output name if needed
  cluster_name <- paste0("C", i)
  run_all_tests(GOdata, cluster_name)
  
  # Optional: Print progress
  cat("Processed cluster", i, "\n")
}



