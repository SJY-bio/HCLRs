library(Seurat)
library(qs)
library(reticulate)
library(hdf5r)
library(sceasy)
library(BiocParallel)
library(data.table)
library(tidyr)
loompy <- reticulate::import('loompy')
# register(MulticoreParam(workers = 4, progressbar = TRUE)) 
IPF <- dir("GSE176078/sc_rds/")
IPF <- sub("\\.rds$", "", IPF)
sp = 1

for (sp in 1:length(IPF)) {
  # Seurat to AnnData
  IPF_seurat = readRDS(paste0("sc_rds/",dir("GSE176078/sc_rds/")[sp]))
  IPF_seurat <- JoinLayers(IPF_seurat)
  IPF_seurat <- FindVariableFeatures(IPF_seurat,verbose = F) %>%
    NormalizeData(verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
  IPF_seurat[["RNA"]] <- as( IPF_seurat[["RNA"]], "Assay")
  sceasy::convertFormat( IPF_seurat, from="seurat", to="anndata", assay = "RNA", main_layer = "count",
                         outFile=paste0("GSE176078/sc_h5ad/",IPF[sp],".h5ad"))
  meta <- IPF_seurat@meta.data
  df_meta <- tibble(
    Cell      = rownames(meta),
    cell_type = meta$cell.type
  )
  fwrite(df_meta,file = paste0("GSE176078/sc_h5ad/",IPF[sp],"_meta.txt"),quote = F,sep = "\t")
}








