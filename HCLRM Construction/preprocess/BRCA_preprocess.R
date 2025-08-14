##### SC #####
library(Seurat)
library(rhdf5)
library(tidyverse)
library(Matrix)
library(dplyr)
library(tidyr)
library(data.table)

# load data
seRNA <- Read10X("D:/ATAC/HCLRs/GSE176078/MTX",gene.column = 1)
seRNA <- CreateSeuratObject(seRNA,min.cells = 3,min.features = 200)
ann <- read.table("D:/ATAC/HCLRs/GSE176078/BRCA_GSE176078_CellMetainfo_table.tsv",header = T,sep = "\t")
seRNA <- subset(seRNA, cells = ann$Cell)
rownames(ann) <- ann$Cell
ann <- ann[rownames(seRNA@meta.data),]

# sum(rownames(seRNA@meta.data)==ann$Cell)

seRNA@meta.data$Celltype_malignancy <- ann$Celltype..malignancy.
seRNA@meta.data$cell.type <- ann$Celltype..minor.lineage.
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type %in% c("B","Plasma"))] <- "B_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type %in% c("CD4Tn","CD8Tex","Tprolif"))] <- "T_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "M1")] <- "macrophage"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Malignant")] <- "cancer_cell"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Monocyte")] <- "myeloid_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "pDC")] <- "DC_cells"
table(seRNA@meta.data$cell.type)
saveRDS(seRNA,"D:/ATAC/HCLRs/GSE176078/GSE176078.rds")

seRNA@meta.data$Sample <- ann$Patient
seRNA <- subset(seRNA,Celltype_malignancy %in% c("Immune cells","Malignant cells"))

sampledata <- unique(seRNA$Sample)
for(x in sampledata){
  sub <- subset(seRNA, Sample == x)
  print(x)
  print(table(sub@meta.data$Celltype_malignancy))
  saveRDS(sub,paste("D:/ATAC/HCLRs/GSE176078/sc_rds/", x, ".rds", sep = ""))
}
### test
seRNA <- readRDS("D:\\ATAC\\HCLRs\\GSE176078\\sc_rds\\CID3586.rds")
seRNA <- CreateSeuratObject(seRNA,min.cells = 3,min.features = 200)



##### ST #####
root_dir <- "D:/ATAC/HCLRs/GSE176078/spatial"

# get the first level subdirectories (each sample folder)
sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
sample_path = sample_dirs[3]
for (sample_path in sample_dirs) {
  # sample name
  sample_name <- basename(sample_path)
  
  # the folder of the text matrix (note the name should be consistent with the actual one)
  mat_dir <- file.path(sample_path, paste0(sample_name, "_filtered_count_matrix"))
  if (!dir.exists(mat_dir)) {
    warning("目录不存在：", mat_dir)
    next
  }
  
  # read the 10X text matrix
  message("Reading 10X matrix for sample ", sample_name, " ...")
  mat <- Read10X(data.dir = mat_dir,gene.column = 1)

  # output H5 file path: note the double underscore
  out_h5 <- file.path(sample_path, paste0(sample_name, "__filtered_feature_bc_matrix.h5"))
  
  # write as 10X HDF5 format
  message("Writing H5 for sample ", sample_name, " to ", out_h5, " ...")
  write10xCounts(path = out_h5,x = mat)
}
### test
# seRNA <- Read10X_h5("D:\\ATAC\\HCLRs\\GSE176078\\spatial\\1142243F\\1142243F__filtered_feature_bc_matrix.h5")
# seRNA <- CreateSeuratObject(seRNA,min.cells = 3,min.features = 200)



#################### GSE176078 空转细胞类型注释 #######################

library(spacexr)
library(Seurat)
library(dplyr)
library(Matrix)

# set the path
spatial_root <- "D:/ATAC/HCLRs/GSE176078/spatial"  # the spatial data directory
output_dir <- "D:/ATAC/HCLRs/GSE176078/spatial_RCTD"  # the output directory

# 创建输出目录
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# load the single cell reference data
sc_ref <- readRDS("D:/ATAC/HCLRs/GSE176078/GSE176078.rds")

# extract the single cell expression matrix and cell type annotation
sc_counts <- GetAssayData(sc_ref, slot = "counts")
sc_meta <- sc_ref@meta.data

# check the cell type column name (adjust according to the actual data)
celltype_col <- "cell.type"  # modify to the actual cell type column name
if (!celltype_col %in% colnames(sc_meta)) {
  stop("Cell type column not found in single cell metadata")
}

# filter the low-abundance cell types (at least 25 cells)
cell_counts <- table(sc_meta[[celltype_col]])
keep_celltypes <- names(cell_counts)[cell_counts >= 25]
sc_meta <- sc_meta[sc_meta[[celltype_col]] %in% keep_celltypes, ]
sc_counts <- sc_counts[, rownames(sc_meta)]

# prepare the Reference object
cluster <- as.factor(sc_meta[[celltype_col]])
names(cluster) <- colnames(sc_counts)
nUMI <- colSums(sc_counts)
reference <- Reference(as.matrix(sc_counts), cluster, nUMI, min_UMI = 50)

# get all the spatial sample directories
spatial_dirs <- list.dirs(spatial_root, recursive = FALSE)

# process each spatial sample
spatial_dir = spatial_dirs[1]
for (spatial_dir in spatial_dirs) {
  sample_id <- basename(spatial_dir)
  cat("Processing sample:", sample_id, "\n")
  
  # build the file path
  h5_file <- file.path(spatial_dir, paste0(sample_id, "__filtered_feature_bc_matrix.h5"))
  positions_file <- file.path(spatial_dir, "spatial/tissue_positions_list.csv")
  json_file <- file.path(spatial_dir, "spatial/scalefactors_json.json")
  
  # read the spatial data
  st_counts <- Read10X_h5(h5_file)
  st_coords <- read.csv(positions_file, header = FALSE, 
                        col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
  
  # filter the tissue area
  st_coords <- st_coords[st_coords$tissue == 1, ]
  
  # ensure the coordinates match the expression matrix
  common_barcodes <- intersect(colnames(st_counts), st_coords$barcode)
  st_counts <- st_counts[, common_barcodes]
  st_coords <- st_coords[match(common_barcodes, st_coords$barcode), ]
  
  # create the spatial coordinate matrix
  coord_mat <- st_coords[, c("imagerow", "imagecol")]
  rownames(coord_mat) <- st_coords$barcode
  colnames(coord_mat) <- c("x", "y")
  
  # create the SpatialRNA object
  nUMI_st <- colSums(st_counts)
  spatial_rna <- SpatialRNA(coords = coord_mat, 
                            counts = st_counts, 
                            nUMI = nUMI_st)
  
  # run RCTD
  RCTD <- create.RCTD(spatial_rna, reference, max_cores = 8)
  RCTD_results <- run.RCTD(RCTD, doublet_mode = "full")
  
  # extract the results
  weights <- RCTD_results@results$weights
  norm_weights <- normalize_weights(weights)  # normalize the weights
  
  # get the main cell types
  cell_types <- colnames(norm_weights)
  pred_celltypes <- apply(norm_weights, 1, function(x) {
    cell_types[which.max(x)]
  })
  
  # create the output data frame
  result_df <- data.frame(
    barcode = names(pred_celltypes),
    celltype = pred_celltypes,
    stringsAsFactors = FALSE
  )
  
  # save the results
  output_file <- file.path(output_dir, paste0(sample_id, "_RCTD_celltypes.txt"))
  write.table(result_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Saved cell type annotations for", nrow(result_df), "spots to", output_file, "\n")
  cat("Completed processing for sample:", sample_id, "\n\n")
}






