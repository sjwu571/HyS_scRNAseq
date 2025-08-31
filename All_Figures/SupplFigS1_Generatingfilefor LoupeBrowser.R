install.packages("remotes")
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()

library(Seurat)
library(loupeR)

# Reorder factor levels so Loupe shows them in your preferred order
seurat_integrated$sample <- factor(
  seurat_integrated$sample,
  levels = c("Methanol_New", "ACME", "Live_cells")
)
# Optional: reorder columns if desired
seurat_integrated <- seurat_integrated[, order(seurat_integrated$sample)]

# Minimal call to generate the H5 (it may throw an error at the end, ignore it)
create_loupe_from_seurat(seurat_integrated) ##If you just start here I think it just leaves the column order as default

# Locate the newest H5 file in the temp folder
tmp_dir <- tempdir()
h5_files <- list.files(tmp_dir, pattern = "\\.h5$", full.names = TRUE)
tmp_h5 <- h5_files[which.max(file.info(h5_files)$ctime)]

cat("Temporary H5 created at:\n", tmp_h5, "\n")
##Copy the path printed — that’s the H5 file Louper needs.

##Open Terminal and run:
"/Users/danielledejong/Library/Application Support/org.R-project.R/R/loupeR/louper" create \
--input="/full/path/to/the/copied/fileXXXX.h5" \ ##Replace /full/path/to/the/copied/fileXXXX.h5 with the path from R.
--output="my_atlas_ordered.cloupe" #or whatever you want your file name to be

##The .cloupe will be created in your current folder (or specify a full path if you prefer).

"/Users/danielledejong/Library/Application Support/org.R-project.R/R/loupeR/louper" create \
--input="/var/folders/b2/ngxm1ykj0yz3fshy20jy0hsr0000gn/T//Rtmp6c4XdE/fileb7dd1f6dac5b.h5" \
--output="my_atlas_ordered_v2.cloupe" #or whatever you want your file name to be
