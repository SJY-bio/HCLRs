########### LUAD_preprocess(GSE127465+GSE179572) ###########
##### SC #####
library(Seurat)
library(rhdf5)
library(tidyverse)
library(Matrix)
library(dplyr)
library(tidyr)
library(data.table)

# load packages
library(R.utils)
library(data.table)
library(Matrix)
library(Seurat)
library(stringr)

# create list to store Seurat objects
seurat_list <- list()

# process each compressed file
setwd("D:\\ATAC\\HCLRs\\GSE127465\\unzipped")
new_file = dir("D:\\ATAC\\HCLRs\\GSE127465\\unzipped")[1]
for (new_file in dir("D:\\ATAC\\HCLRs\\GSE127465\\unzipped")) {
  sample_id <- strsplit(new_file,".tsv")[[1]]
  # read TSV file
  dt <- fread(new_file, header = TRUE)
  
  # convert to sparse matrix (rows=cells, columns=genes)
  counts_matrix <- as(as.matrix(dt[, -1]), "sparseMatrix")
  rownames(counts_matrix) <- dt$barcode
  counts_matrix <- t(counts_matrix)  # transpose to Seurat required format (rows=genes, columns=cells)
  
  # create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts_matrix,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  # add to list
  seurat_list[[sample_id]] <- seurat_obj
  
  # clean up temporary variables
  rm(dt, counts_matrix)
  gc()
}

# merge all Seurat objects
seRNA <- merge(seurat_list[[1]], 
                       y = unlist(seurat_list[-1]), 
                       add.cell.ids = names(seurat_list))
seRNA <- JoinLayers(seRNA)
ann <- fread("D:/ATAC/HCLRs/GSE127465/NSCLC_GSE127465_CellMetainfo_table.tsv",header = T,check.names = T)
ann <- na.omit(ann)
ann <- as.data.frame(ann)
seRNA <- subset(seRNA, cells = ann$Cell)
rownames(ann) <- ann$Cell
ann <- ann[rownames(seRNA@meta.data),]

# sum(rownames(seRNA@meta.data)==ann$Cell)

seRNA@meta.data$Celltype_malignancy <- ann$Celltype..malignancy.
seRNA@meta.data$cell.type <- ann$Celltype..minor.lineage.
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type %in% c("B","Plasma"))] <- "B_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type %in% c("CD4Tn","CD8Tex","Th2"))] <- "T_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "M2")] <- "macrophage"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Malignant")] <- "cancer_cell"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Monocyte")] <- "myeloid_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "cDC2")] <- "DC_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "NK")] <- "NK_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Mast")] <- "mast_cells"
seRNA@meta.data$cell.type[which(seRNA@meta.data$cell.type == "Neutrophils")] <- "neutrophils"

table(seRNA@meta.data$cell.type)
saveRDS(seRNA,"D:/ATAC/HCLRs/GSE127465/GSE127465_sc.rds")

seRNA@meta.data$Sample <- ann$Patient
seRNA <- subset(seRNA,Celltype_malignancy %in% c("Immune cells","Malignant cells"))

sampledata <- unique(seRNA$Sample)
for(x in sampledata){
  sub <- subset(seRNA, Sample == x)
  print(x)
  print(table(sub@meta.data$Celltype_malignancy))
  print(table(sub@meta.data$cell.type))
  saveRDS(sub,paste("D:/ATAC/HCLRs/GSE127465/sc_rds/", x, ".rds", sep = ""))
}
### test
# seRNA <- readRDS("D:\\ATAC\\HCLRs\\GSE117570_LUAD\\sc_rds\\P1.rds")


##### ST #####
library(spacexr)
library(Seurat)
library(dplyr)
library(Matrix)

# set path
spatial_root <- "D:/ATAC/HCLRs/GSE179572/spatial"  # spatial data directory
output_dir <- "D:/ATAC/HCLRs/GSE179572/spatial_RCTD"  # output directory

# create output directory
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# load single cell reference data
sc_ref <- readRDS("D:/ATAC/HCLRs/GSE127465/GSE127465_sc.rds")

# extract single cell expression matrix and cell type annotation
sc_counts <- GetAssayData(sc_ref, slot = "counts")
sc_meta <- sc_ref@meta.data

# check cell type column name (adjust according to actual data)
celltype_col <- "cell.type"  # modify to actual cell type column name
if (!celltype_col %in% colnames(sc_meta)) {
  stop("Cell type column not found in single cell metadata")
}

# filter low abundance cell types (at least 25 cells)
cell_counts <- table(sc_meta[[celltype_col]])
keep_celltypes <- names(cell_counts)[cell_counts >= 25]
sc_meta <- sc_meta[sc_meta[[celltype_col]] %in% keep_celltypes, ]
sc_counts <- sc_counts[, rownames(sc_meta)]

# prepare Reference object
cluster <- as.factor(sc_meta[[celltype_col]])
names(cluster) <- colnames(sc_counts)
nUMI <- colSums(sc_counts)
reference <- Reference(as.matrix(sc_counts), cluster, nUMI, min_UMI = 50)

# get all spatial sample directories
spatial_dirs <- list.dirs(spatial_root, recursive = FALSE)

# process each spatial sample
spatial_dir = spatial_dirs[1]
for (spatial_dir in spatial_dirs) {
  sample_id <- basename(spatial_dir)
  cat("Processing sample:", sample_id, "\n")
  
  # build file path
  h5_file <- file.path(spatial_dir, paste0(sample_id, "__filtered_feature_bc_matrix.h5"))
  positions_file <- file.path(spatial_dir, "spatial/tissue_positions_list.csv")
  json_file <- file.path(spatial_dir, "spatial/scalefactors_json.json")
  
  # read spatial data
  st_counts <- Read10X_h5(h5_file)
  st_coords <- read.csv(positions_file, header = FALSE, 
                        col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
  
  # filter tissue area
  st_coords <- st_coords[st_coords$tissue == 1, ]
  
  # ensure coordinates match expression matrix
  common_barcodes <- intersect(colnames(st_counts), st_coords$barcode)
  st_counts <- st_counts[, common_barcodes]
  st_coords <- st_coords[match(common_barcodes, st_coords$barcode), ]
  
  # create spatial coordinate matrix
  coord_mat <- st_coords[, c("imagerow", "imagecol")]
  rownames(coord_mat) <- st_coords$barcode
  colnames(coord_mat) <- c("x", "y")
  
  # create SpatialRNA object
  nUMI_st <- colSums(st_counts)
  spatial_rna <- SpatialRNA(coords = coord_mat, 
                            counts = st_counts, 
                            nUMI = nUMI_st)
  
  # run RCTD
  RCTD <- create.RCTD(spatial_rna, reference, max_cores = 8)
  RCTD_results <- run.RCTD(RCTD, doublet_mode = "full")
  
  # extract results
  weights <- RCTD_results@results$weights
  norm_weights <- normalize_weights(weights)  # normalize weights
  
  # get main cell types
  cell_types <- colnames(norm_weights)
  pred_celltypes <- apply(norm_weights, 1, function(x) {
    cell_types[which.max(x)]
  })
  
  # create output data frame
  result_df <- data.frame(
    barcode = names(pred_celltypes),
    celltype = pred_celltypes,
    stringsAsFactors = FALSE
  )
  
  # save results
  output_file <- file.path(output_dir, paste0(sample_id, "_RCTD_celltypes.txt"))
  write.table(result_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Saved cell type annotations for", nrow(result_df), "spots to", output_file, "\n")
  cat("Completed processing for sample:", sample_id, "\n\n")
}







