
# -----------------------------
# 1. Load Required Libraries
# -----------------------------
library(ggplot2)
library(ggridges)
library(readxl)
library(Seurat)
library(CellChat)
library(SpaTalk)

# -----------------------------
# 2. Set Paths (Use Relative Paths for GitHub)
# -----------------------------
data_dir <- "data/raw/"
results_dir <- "results/interaction/"
figures_dir <- "results/figures/"

# Create output directories
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3. Part A: Analyze Near/Far Interactions Using CellChat
# -----------------------------

# Load Seurat object and metadata
setwd(data_dir)
load("SCC_scRNA_exp_scale_seurat_obj.Rdata")
spot_near_far_meta <- read.table("spot_near_far_meta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SCC_near_far_exp <- read.table("st_exp_symbol_near_far.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create Seurat object
SCC_scRNA_naer_far <- CreateSeuratObject(counts = SCC_near_far_exp, min.cells = 3, min.features = 10)
SCC_scRNA_naer_far <- NormalizeData(SCC_scRNA_naer_far)
exp <- as.matrix(SCC_scRNA_naer_far@assays$RNA@data)
spot_near_far_meta$celltype <- gsub("_", " ", spot_near_far_meta$celltype)

# Split expression by Near/Far + Normal
exp_far <- exp[, which(spot_near_far_meta$near_far_mali %in% c("Far", "Normal_Epithelial"))]
meta_far <- spot_near_far_meta[rownames(exp_far), ]
exp_near <- exp[, which(spot_near_far_meta$near_far_mali %in% c("Near", "Normal_Epithelial"))]
meta_near <- spot_near_far_meta[rownames(exp_near), ]

# Load interaction database
inter_all <- read.table("union.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
inter_all$interaction_name <- paste(inter_all$source_genesymbol, inter_all$target_genesymbol, sep = "_")
inter_all <- inter_all[!duplicated(inter_all$interaction_name), c("interaction_name", "ligand" = "source_genesymbol", "receptor" = "target_genesymbol")]

# Build custom CellChat DB
interaction_union <- data.frame(
  interaction_name = inter_all$interaction_name,
  pathway_name = "NULL", ligand = inter_all$ligand, receptor = inter_all$receptor,
  agonist = "NULL", antagonist = "NULL", co_A_receptor = "NULL",
  co_I_receptor = "NULL", interaction_name_2 = "NULL",
  annotation = "NULL", evidence = "NULL"
)
rownames(interaction_union) <- interaction_union$interaction_name

resource_union_DB <- list(
  interaction = interaction_union,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  geneInfo = CellChatDB.human$geneInfo
)

# Initialize result storage
resource22_Near_LRscore1 <- data.frame()

# Run CellChat on Near condition
for (n in 1:21) {
  cat("Processing iteration:", n, "\n")
  
  cellchat_scc <- createCellChat(exp_near, meta_near, group.by = "near_far_mali")
  cellchat_scc@DB <- resource_union_DB
  cellchat_scc <- subsetData(cellchat_scc)
  cellchat_scc <- identifyOverExpressedGenes(cellchat_scc)
  cellchat_scc <- identifyOverExpressedInteractions(cellchat_scc)
  cellchat_scc <- projectData(cellchat_scc, PPI.human)
  cellchat_scc <- computeCommunProb(cellchat_scc, type = "truncatedMean", trim = 0.1, raw.use = TRUE)
  cellchat_scc <- filterCommunication(cellchat_scc, min.cells = 1)
  
  scc_net <- subsetCommunication(cellchat_scc)
  scc_net$LR <- paste(scc_net$ligand, scc_net$receptor, sep = "_")
  scc_net$CC <- paste(scc_net$source, scc_net$target, sep = "_")
  
  # Only keep cross-type interactions
  scc_net_filtered <- scc_net[!(scc_net$CC %in% c("Near_Near", "Normal_Epithelial_Normal_Epithelial")), ]
  if (nrow(scc_net_filtered) == 0) next
  
  # Compute custom interaction score: sum(ligand_exp) * sum(receptor_exp)
  scc_net_filtered$LR_near_score <- sapply(1:nrow(scc_net_filtered), function(i) {
    ligand <- scc_net_filtered$ligand[i]
    receptor <- scc_net_filtered$receptor[i]
    src <- scc_net_filtered$source[i]
    tgt <- scc_net_filtered$target[i]
    
    lig_exp <- exp[rownames(exp) == ligand, colnames(exp) %in% meta_near$cell_ID[meta_near$near_far_mali == src]]
    rec_exp <- exp[rownames(exp) == receptor, colnames(exp) %in% meta_near$cell_ID[meta_near$near_far_mali == tgt]]
    
    if (length(lig_exp) == 0 || length(rec_exp) == 0) return(0)
    return(sum(lig_exp) * sum(rec_exp))
  })
  
  # Add small noise to zero values for log transformation
  scc_net_filtered$log_score <- ifelse(
    scc_net_filtered$LR_near_score == 0,
    log(runif(sum(scc_net_filtered$LR_near_score == 0), 0, 0.001)),
    log(scc_net_filtered$LR_near_score)
  )
  
  scc_net_filtered <- scc_net_filtered[, c(1:6, 12:15)]
  scc_net_filtered$resource <- "Union"
  resource22_Near_LRscore1 <- rbind(resource22_Near_LRscore1, scc_net_filtered)
}

# Save results
write.table(resource22_Near_LRscore1, "union_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 4. Part B: Visualization (Ridge Plots, Boxplots)
# -----------------------------

# Read results with custom scores
resource22_Near_LRscore1 <- read_xlsx("resource23_near_norm_cellchat_LRscore.xlsx")
resource22_Near_LRscore1$log_score <- as.numeric(resource22_Near_LRscore1$log_score)
resource22_Near_LRscore1$prob <- as.numeric(resource22_Near_LRscore1$prob)

# Load color scheme
color22 <- read.table("21resource_color.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
colors <- setNames(color22$color2, color22$resource)

# Ridge plot: log score
pdf(file.path(figures_dir, "resource23_score_new.pdf"), height = 6, width = 4.5)
ggplot(resource22_Near_LRscore1, aes(y = resource, x = log_score, fill = resource)) +
  geom_density_ridges_gradient(scale = 0.6, alpha = 0.7) +
  scale_fill_manual(values = colors) +
  scale_y_discrete(limits = rev(color22$resource)) +
  xlab("CCLR score") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# Ridge plot: probability
pdf(file.path(figures_dir, "resource23_score_new_prob.pdf"), height = 6, width = 4.5)
ggplot(resource22_Near_LRscore1, aes(y = resource, x = prob, fill = resource)) +
  geom_density_ridges_gradient(scale = 0.6, alpha = 0.7) +
  scale_fill_manual(values = colors) +
  scale_y_discrete(limits = rev(color22$resource)) +
  xlab("prob") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# Bar plot: LR count per resource
resource22_sig_LRcount <- data.frame(table(resource22_Near_LRscore1$resource))
pdf(file.path(figures_dir, "resource22_LRcount_bang_new.pdf"), height = 6, width = 3)
ggplot(resource22_sig_LRcount, aes(x = Freq, y = Var1)) +
  geom_col(fill = "#1170AA", width = 0.2) +
  geom_point(color = "#FC7D0B", size = 4) +
  ylab("") + xlab("LR number") +
  theme_bw() + theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(color22$resource)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# -----------------------------
# 5. Part C: Compare HCLRs vs Union (Boxplots)
# -----------------------------

data_all <- read_xlsx("resource23_near_norm_cellchat_LRscore.xlsx")
data_all <- data_all[data_all$resource %in% c("HCLRs", "Union"), ]
data_all$log_score <- as.numeric(data_all$log_score)
data_all$prob <- as.numeric(data_all$prob)

# Boxplot: log score
p_score <- ggplot(data_all, aes(y = log_score, x = resource, fill = resource, color = resource)) +
  scale_fill_manual(values = c("#B52027", "#BAB0AC")) +
  scale_color_manual(values = c("#B52027", "#BAB0AC")) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  ylim(0, 12) +
  theme_classic() +
  xlab("") +
  theme(axis.ticks.length = unit(0.15, "cm"), text = element_text(size = 20)) +
  stat_compare_means(comparisons = list(c("HCLRs", "Union")), label = "p.signif", method = "wilcox.test")

ggsave(p_score, file = file.path(figures_dir, "HCLRs_Union_score.pdf"), width = 6, height = 5)

# Boxplot: probability
p_prob <- ggplot(data_all, aes(y = prob, x = resource, fill = resource, color = resource)) +
  scale_fill_manual(values = c("#B52027", "#BAB0AC")) +
  scale_color_manual(values = c("#B52027", "#BAB0AC")) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +
  ylim(-0.1, 0.3) +
  theme_classic() +
  xlab("") +
  theme(axis.ticks.length = unit(0.15, "cm"), text = element_text(size = 20)) +
  stat_compare_means(comparisons = list(c("HCLRs", "Union")), label = "p.signif", method = "wilcox.test")

ggsave(p_prob, file = file.path(figures_dir, "HCLRs_Union_prob.pdf"), width = 6, height = 5)

# -----------------------------
# 6. Part D: FOV-Level Analysis Using SpaTalk
# -----------------------------

# Load FOV data
cell_meta_lung13 <- read.table(file.path(data_dir, "basic_Lung13_data.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Lung13_xy_meta <- read.table(file.path(data_dir, "Lung13_xy_meta.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Lung13_exp_final <- read.table(file.path(data_dir, "Lung13_exp.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
CC_dist_tumor_other <- read.table(file.path(data_dir, "CC_dist_tumor_other_meta_near_far.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

high_conf_LR <- read.table(file.path(data_dir, "xgboost_confidence_LR_01_8777.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
High_LR_pairs <- data.frame(ligand = high_conf_LR$source, receptor = high_conf_LR$target, species = "Human")

# Loop over 20 FOVs
for (i in 1:20) {
  message("Processing FOV: ", i)
  
  meta_20 <- Lung13_xy_meta[Lung13_xy_meta$fov == i, ]
  dist_fov <- CC_dist_tumor_other[CC_dist_tumor_other$to %in% meta_20$cell_ID & CC_dist_tumor_other$from %in% meta_20$cell_ID, ]
  dist_fov <- dist_fov[order(dist_fov$distance), ]
  
  n <- nrow(dist_fov)
  dist_fov$near_far_other <- "Other"
  dist_fov$near_far_other[1:ceiling(n/4)] <- "Near"
  dist_fov$near_far_other[floor(3*n/4):n] <- "Far"
  
  # Extract Near/Far cells with cell types
  extract_cells <- function(df, group) {
    from <- df[df$near_far_other == group, c("from", "from_ct")]
    to <- df[df$near_far_other == group, c("to", "to_ct")]
    colnames(from) <- colnames(to) <- c("cell", "celltype")
    unique(rbind(from, to))
  }
  
  near_cells <- extract_cells(dist_fov, "Near")
  far_cells  <- extract_cells(dist_fov, "Far")
  
  # Add coordinates
  near_with_xy <- merge(near_cells, meta_20, by.x = "cell", by.y = "cell_ID")
  far_with_xy  <- merge(far_cells,  meta_20, by.x = "cell", by.y = "cell_ID")
  
  # Subset expression
  exp_near <- Lung13_exp_final[, colnames(Lung13_exp_final) %in% near_with_xy$cell]
  exp_far  <- Lung13_exp_final[, colnames(Lung13_exp_final) %in% far_with_xy$cell]
  
  # Create SpaTalk objects
  meta_near <- data.frame(cell = near_with_xy$cell, x = near_with_xy$sdimx, y = near_with_xy$sdimy)
  meta_far  <- data.frame(cell = far_with_xy$cell,  x = far_with_xy$sdimx,  y = far_with_xy$sdimy)
  
  obj_near <- createSpaTalk(exp_near, meta_near, if_st_is_sc = TRUE, spot_max_cell = 1, celltype = near_with_xy$celltype, species = "Human")
  obj_far  <- createSpaTalk(exp_far,  meta_far,  if_st_is_sc = TRUE, spot_max_cell = 1, celltype = far_with_xy$celltype, species = "Human")
  
  obj_near <- find_lr_path(obj_near, lrpairs = High_LR_pairs, pathways = pathways)
  obj_far  <- find_lr_path(obj_far,  lrpairs = High_LR_pairs, pathways = pathways)
  
  # Define tumor-other interactions
  ct_types_near <- setdiff(unique(obj_near@meta[["rawmeta"]][["celltype"]]), "tumor")
  ct_types_far  <- setdiff(unique(obj_far@meta[["rawmeta"]][["celltype"]]),  "tumor")
  
  for (ct in ct_types_near) {
    # Tumor -> CT
    res1 <- tryCatch({ dec_cci(obj_near, "tumor", ct, n_neighbor = 2)@lrpair }, error = function(e) NULL)
    # CT -> Tumor
    res2 <- tryCatch({ dec_cci(obj_near, ct, "tumor", n_neighbor = 2)@lrpair }, error = function(e) NULL)
    
    combined <- na.omit(rbind(res1, res2))
    if (!is.null(combined) && nrow(combined) > 0) {
      combined$near_far <- "Near"
      combined$fov <- paste0("fov", i)
      combined$CC_inter <- paste0("tumor--", ct)
      combined$resource <- "HCLRs"
      write.table(combined, file.path(results_dir, "Near_fov_union", paste0("fov_", i, "_tumor_", ct, "_near.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  
  for (ct in ct_types_far) {
    res1 <- tryCatch({ dec_cci(obj_far, "tumor", ct, n_neighbor = 2)@lrpair }, error = function(e) NULL)
    res2 <- tryCatch({ dec_cci(obj_far, ct, "tumor", n_neighbor = 2)@lrpair }, error = function(e) NULL)
    
    combined <- na.omit(rbind(res1, res2))
    if (!is.null(combined) && nrow(combined) > 0) {
      combined$near_far <- "Far"
      combined$fov <- paste0("fov", i)
      combined$CC_inter <- paste0("tumor--", ct)
      combined$resource <- "HCLRs"
      write.table(combined, file.path(results_dir, "Far_fov_union", paste0("fov_", i, "_tumor_", ct, "_far.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}

message("✅ All analyses completed and results saved.")


# -----------------------------
# -----------------------------
# 1. Load Required Libraries
# -----------------------------
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(paletteer)
library(ggpubr)

# -----------------------------
# 2. Set Working Directory and Paths
# -----------------------------
# Use relative paths for portability
data_dir <- "data/raw"
results_dir <- "results/interaction"
figures_dir <- "results/figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

setwd(data_dir)

# -----------------------------
# 3. Load and Preprocess Data
# -----------------------------
# Load Seurat object and expression matrix
load("SCC_scRNA_exp_scale_seurat_obj.Rdata")
spot_near_far_meta <- read.table("scRNA_celltype_metadata_filter_tumor_Normal.txt", 
                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SCC_near_far_exp <- read.table("scRNA_data_exp_Tumor_Normal.txt", 
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create Seurat object and normalize
SCC_scRNA_naer_far <- CreateSeuratObject(counts = SCC_near_far_exp, min.cells = 3, min.features = 10)
SCC_scRNA_naer_far <- NormalizeData(SCC_scRNA_naer_far)
exp <- as.matrix(SCC_scRNA_naer_far@assays$RNA@data)

# Fix cell type naming
spot_near_far_meta$celltype <- gsub("_", " ", spot_near_far_meta$celltype)

# -----------------------------
# 4. Define Key Cell Interactions of Interest
# -----------------------------
target_interactions <- c(
  'Endothelial Cell--Fibroblast', 'Fibroblast--Endothelial Cell',
  'Malignant Epithelial--B Cell', 'B Cell--Malignant Epithelial',
  'Malignant Epithelial--DCs', 'DCs--Malignant Epithelial',
  'Malignant Epithelial--Endothelial Cell', 'Endothelial Cell--Malignant Epithelial',
  'Malignant Epithelial--Fibroblast', 'Fibroblast--Malignant Epithelial',
  'Malignant Epithelial--Mono/Macro', 'Mono/Macro--Malignant Epithelial',
  'Normal Epithelial--Mono/Macro', 'Mono/Macro--Normal Epithelial',
  'Normal Epithelial--DCs', 'DCs--Normal Epithelial',
  'Normal Epithelial--Endothelial Cell', 'Endothelial Cell--Normal Epithelial',
  'Normal Epithelial--Fibroblast', 'Fibroblast-Normal Epithelial',
  'Normal Epithelial--Malignant Epithelial', 'Malignant Epithelial--Normal Epithelial'
)

# Define reverse patterns for consistent orientation
reverse_patterns <- c(
  "Fibroblast--Endothelial Cell",
  "B Cell--Malignant Epithelial",
  "DCs--Malignant Epithelial",
  "Endothelial Cell--Malignant Epithelial",
  "Fibroblast--Malignant Epithelial",
  "Mono/Macro--Malignant Epithelial",
  "Normal Epithelial--Mono/Macro",
  "DCs--Normal Epithelial",
  "Endothelial Cell--Normal Epithelial",
  "Fibroblast-Normal Epithelial",
  "Malignant Epithelial--Normal Epithelial"
)

# -----------------------------
# 5. Run CellChat with Custom LR Database
# -----------------------------
# Load high-confidence LR list
inter_all <- read.table('HCLR_new.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
inter_all$interaction_name <- paste(inter_all$source_genesymbol, inter_all$target_genesymbol, sep = "_")
inter_all <- inter_all[!duplicated(inter_all$interaction_name), ]
inter_all <- inter_all[, c("interaction_name", "ligand" = "source_genesymbol", "receptor" = "target_genesymbol")]

# Build custom CellChat DB
interaction_union <- data.frame(
  interaction_name = inter_all$interaction_name,
  pathway_name = "NULL", ligand = inter_all$ligand, receptor = inter_all$receptor,
  agonist = "NULL", antagonist = "NULL", co_A_receptor = "NULL",
  co_I_receptor = "NULL", interaction_name_2 = "NULL",
  annotation = "NULL", evidence = "NULL"
)
rownames(interaction_union) <- interaction_union$interaction_name

resource_union_DB <- list(
  interaction = interaction_union,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  geneInfo = CellChatDB.human$geneInfo
)

# Initialize result storage
scc_net_near1 <- NULL

# Run CellChat using level2 cell types
cellchat_scc <- createCellChat(object = exp, meta = spot_near_far_meta, group.by = "level2_celltype")
cellchat_scc@DB <- resource_union_DB
cellchat_scc <- subsetData(cellchat_scc)
cellchat_scc <- identifyOverExpressedGenes(cellchat_scc)
cellchat_scc <- identifyOverExpressedInteractions(cellchat_scc)
cellchat_scc <- projectData(cellchat_scc, PPI.human)
cellchat_scc <- computeCommunProb(cellchat_scc, type = "truncatedMean", trim = 0.1, raw.use = TRUE)
cellchat_scc <- filterCommunication(cellchat_scc, min.cells = 1)

# Extract communication network
scc_net <- subsetCommunication(cellchat_scc)
scc_net$inter_cell <- paste(scc_net$source, scc_net$target, sep = "--")
scc_net$LR_pair <- paste(scc_net$ligand, scc_net$receptor, sep = "_")
scc_net$resource <- 'HCLRs'

scc_net_near1 <- scc_net

# Save results
write.table(scc_net_near1, "../HCLRs_new_result.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# -----------------------------
# 6. Combine Results from All Resources
# -----------------------------
all_result <- read.table("all_22resource_result_SCC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_result <- all_result[all_result$resource != 'HCLRs', ]
all_result01 <- read.table("HCLRs_new_result.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

all_result <- rbind(all_result, all_result01)
all_result <- all_result[all_result$inter_cell %in% target_interactions, ]

# Reverse direction for consistent comparison
all_result <- all_result %>%
  mutate(inter_cell = ifelse(inter_cell %in% reverse_patterns,
                             gsub("(.*)--(.*)", "\\2--\\1", inter_cell),
                             inter_cell))

# -----------------------------
# 7. Plot 1: Number of LR Pairs (Boxplot)
# -----------------------------
all_result$inter_cell <- gsub("--", " | ", all_result$inter_cell)
all_result$num <- 1
all_11ct_LRnum <- aggregate(num ~ inter_cell + resource, data = all_result, sum)
colnames(all_11ct_LRnum) <- c("cellinter", "resource", "LRnum")

# Load color scheme
color22 <- read.table("21resource_color.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
color_22 <- setNames(color22$color2, color22$resource)
colors_fill <- color_22
colors_point <- paletteer_d("ggsci::planetexpress_futurama")[1:length(unique(all_11ct_LRnum$cellinter))]

# Generate boxplot
pdf(file.path(figures_dir, "LRnum_22resouce_11inter_new.pdf"), width = 13, height = 7)
ggplot(all_11ct_LRnum, aes(x = resource, y = LRnum, fill = resource)) +
  geom_boxplot(alpha = 0.9, color = "#49525E", outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", color = "#49525E") +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_fill) +
  geom_point(aes(color = cellinter), size = 3, shape = 5, stroke = 1, alpha = 0.8) +
  geom_line(aes(group = cellinter), color = 'gray', lwd = 0.1) +
  scale_colour_manual(values = colors_point) +
  theme_bw() +
  ylab("Number of LRs") +
  xlab("") +
  scale_x_discrete(limits = color22$resource) +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15))
dev.off()

# Statistical test
wil11 <- compare_means(LRnum ~ resource, data = all_11ct_LRnum)
write.table(wil11, "new/wilxon_22resouce_11inter.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 8. Plot 2: Total Communication Probability
# -----------------------------
all_result$inter_cell <- gsub(" | ", "./", all_result$inter_cell)
all_11ct_prob <- aggregate(prob ~ inter_cell + resource, data = all_result, sum)
colnames(all_11ct_prob) <- c("cellinter", "resource", "prob")
all_11ct_prob$cellinter <- gsub("./", " | ", all_11ct_prob$cellinter)

pdf(file.path(figures_dir, "prob_22resouce_11inter_new.pdf"), width = 13, height = 7)
ggplot(all_11ct_prob, aes(x = resource, y = prob, fill = resource)) +
  geom_boxplot(alpha = 0.9, color = "#49525E", outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", color = "#49525E") +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_fill) +
  geom_point(aes(color = cellinter), size = 3, shape = 5, stroke = 1, alpha = 0.8) +
  geom_line(aes(group = cellinter), color = 'gray', lwd = 0.1) +
  scale_colour_manual(values = colors_point) +
  theme_bw() +
  ylab("Prob") +
  xlab("") +
  scale_x_discrete(limits = color22$resource) +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15))
dev.off()

# Statistical test
wil11_prob <- compare_means(prob ~ resource, data = all_11ct_prob)
write.table(wil11_prob, "new_prob/wilxon_22resouce_11inter.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 9. Visualization: Ridge Plots for HCLRs vs Union
# -----------------------------
# Load full results with log scores
resource22_Near_LRscore1 <- read_xlsx("../resource23_near_norm_cellchat_LRscore.xlsx")
resource22_Near_LRscore1$log_score <- as.numeric(resource22_Near_LRscore1$log_score)
resource22_Near_LRscore1$prob <- as.numeric(resource22_Near_LRscore1$prob)

# Load color scheme
colors <- setNames(color22$color2, color22$resource)

# Ridge plot: log score
pdf(file.path(figures_dir, "resource23_score_new.pdf"), height = 6, width = 4.5)
ggplot(resource22_Near_LRscore1, aes(y = resource, x = log_score, fill = resource)) +
  geom_density_ridges_gradient(scale = 0.6, alpha = 0.7) +
  scale_fill_manual(values = colors) +
  scale_y_discrete(limits = rev(color22$resource)) +
  xlab("CCLR score") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# Ridge plot: probability
pdf(file.path(figures_dir, "resource23_score_new_prob.pdf"), height = 6, width = 4.5)
ggplot(resource22_Near_LRscore1, aes(y = resource, x = prob, fill = resource)) +
  geom_density_ridges_gradient(scale = 0.6, alpha = 0.7) +
  scale_fill_manual(values = colors) +
  scale_y_discrete(limits = rev(color22$resource)) +
  xlab("prob") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# Bar plot: number of significant LRs per resource
lr_count <- data.frame(table(resource22_Near_LRscore1$resource))
pdf(file.path(figures_dir, "resource22_LRcount_bang_new.pdf"), height = 6, width = 3)
ggplot(lr_count, aes(x = Freq, y = Var1)) +
  geom_col(fill = "#1170AA", width = 0.2) +
  geom_point(color = "#FC7D0B", size = 4) +
  ylab("") + xlab("LR number") +
  theme_bw() + theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(color22$resource)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
dev.off()

# -----------------------------
# 10. Final Comparison: HCLRs vs Union
# -----------------------------
data_all <- read_xlsx("../resource23_near_norm_cellchat_LRscore.xlsx")
data_all <- subset(data_all, resource %in% c("HCLRs", "Union"))
data_all$log_score <- as.numeric(data_all$log_score)
data_all$prob <- as.numeric(data_all$prob)

# Boxplot: log score
p_score <- ggplot(data_all, aes(y = log_score, x = resource, fill = resource, color = resource)) +
  scale_fill_manual(values = c("#B52027", "#BAB0AC")) +
  scale_color_manual(values = c("#B52027", "#BAB0AC")) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  ylim(0, 12) +
  theme_classic() +
  xlab("") +
  theme(axis.ticks.length = unit(0.15, "cm"), text = element_text(size = 20)) +
  stat_compare_means(comparisons = list(c("HCLRs", "Union")), label = "p.signif", method = "wilcox.test")

ggsave(p_score, file = file.path(figures_dir, "HCLRs_Union_score.pdf"), width = 6, height = 5)

# Boxplot: probability
p_prob <- ggplot(data_all, aes(y = prob, x = resource, fill = resource, color = resource)) +
  scale_fill_manual(values = c("#B52027", "#BAB0AC")) +
  scale_color_manual(values = c("#B52027", "#BAB0AC")) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +
  ylim(-0.1, 0.3) +
  theme_classic() +
  xlab("") +
  theme(axis.ticks.length = unit(0.15, "cm"), text = element_text(size = 20)) +
  stat_compare_means(comparisons = list(c("HCLRs", "Union")), label = "p.signif", method = "wilcox.test")

ggsave(p_prob, file = file.path(figures_dir, "HCLRs_Union_prob.pdf"), width = 6, height = 5)

message("✅ All analyses completed. Results saved in 'results/' directory.")