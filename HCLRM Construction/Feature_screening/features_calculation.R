##### OmniPath #####
setwd(dir = "D:/ATAC/HCLRs/OmniPath")
library(data.table)
LR <- data.frame()
for (i in 1:length(dir())) {
  a <- fread(dir()[i],header = T,sep = "\t")
  LR <- rbind(LR,a[,c(1:2)])
}
LR <- unique(LR)
Ligand <- data.frame(Ligand = unique(LR$source))
Recepter <- data.frame(Recepter = unique(LR$target))
fwrite(Ligand,file = "D:/ATAC/HCLRs/Ligand_unique.txt",quote = F,sep = "\t",col.names = F)
fwrite(Recepter,file = "D:/ATAC/HCLRs/Recepter_unique.txt",quote = F,sep = "\t",col.names = F)




xgboost <- fread("D:\\ATAC\\HCLRs\\xgboost_confidence_LR_01_8777.txt",header = T,sep = "\t")

##### CellChatDB_conservatism_3species #####
load("D:/ATAC/HCLRs/Species/CellChatDB.human.rda")
load("D:/ATAC/HCLRs/Species/CellChatDB.mouse.rda")
load("D:/ATAC/HCLRs/Species/CellChatDB.zebrafish.rda")
human <- CellChatDB.human$interaction
mouse <- CellChatDB.mouse$interaction
zebrafish <- CellChatDB.zebrafish$interaction
inter_species <- human[which(human$interaction_name%in%intersect(intersect(human$interaction_name,mouse$interaction_name),zebrafish$interaction_name)),]
rm(CellChatDB.human,CellChatDB.mouse,CellChatDB.zebrafish,human,mouse,zebrafish)

inter_species$interaction_name <- paste0(inter_species$ligand,"--",inter_species$receptor)
inter_species <- inter_species$interaction_name

xgboost$Interspecies_conservatism <- 0
xgboost[which(xgboost$pair %in% inter_species),"Interspecies_conservatism"] <- 1

##### KEGG_TF_class #####
KEGG <- fread("D:/ATAC/HCLRs/all_axis_filter_out.csv",header = T,sep = ",")[,1:3]
KEGG$pair <- paste0(KEGG$Ligand_Symbol,"--",KEGG$Receptor_Symbol)
KEGG <- unique(KEGG[,c("pair","TF_Symbol")])
library(dplyr)
KEGG_simplified <- KEGG %>%
  group_by(pair) %>%
  summarise(
    TF_list = paste(unique(TF_Symbol), collapse = ";"),  # 把所有下游 TF 用分号串起来
    KEGG_TF_count = n_distinct(TF_Symbol)                     # 计算每对有多少个不同的TF
  ) %>% ungroup()
KEGG_simplified$KEGG_TF_class <- 1
xgboost <- merge(xgboost,KEGG_simplified[,c(1,3,4)],by = "pair",all.x = T)
xgboost[is.na(xgboost$KEGG_TF_class),"KEGG_TF_class"] <- 0
xgboost[is.na(xgboost$KEGG_TF_count),"KEGG_TF_count"] <- 0

cellcall <- fread("D:/ATAC/HCLRs/new_ligand_receptor_TFs.txt",header = T,sep = "\t")[,6:8]
cellcall$pair <- paste0(cellcall$Ligand_Symbol,"--",cellcall$Receptor_Symbol)
cellcall <- unique(cellcall[,c("pair","TF_Symbol")])
library(dplyr)
cellcall_simplified <- cellcall %>%
  group_by(pair) %>%
  summarise(
    TF_list = paste(unique(TF_Symbol), collapse = ";"),  # 把所有下游 TF 用分号串起来
    cellcall_TF_count = n_distinct(TF_Symbol)                     # 计算每对有多少个不同的TF
  ) %>% ungroup()
cellcall_simplified$cellcall_TF_class <- 1
xgboost <- merge(xgboost,cellcall_simplified[,c(1,3,4)],by = "pair",all.x = T)
xgboost[is.na(xgboost$cellcall_TF_class),"cellcall_TF_class"] <- 0
xgboost[is.na(xgboost$cellcall_TF_count),"cellcall_TF_count"] <- 0


##### Cor_GTEx #####
library(cmapR)
my_ds = parse_gctx("D:/ATAC/HCLRs/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct")
exp <- my_ds@mat

# 转换行名
library(biomaRt)
# a) 去掉行名里的版本号
ens_ids <- rownames(exp)
ens_stripped <- sub("\\..*$", "", ens_ids)  
# b) 从 Ensembl BioMart 拉基因名
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(ens_stripped),
  mart       = mart
)
# 去掉没有对应 symbol 的
mapping <- mapping[mapping$hgnc_symbol != "", ]
# c) 构建一个从原行名 → symbol 的向量
#    注意如果一个 ENSEMBL 对应多个 symbol，这里只取第一个
mapping <- mapping[!duplicated(mapping$ensembl_gene_id), ]
emap <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
# d) 筛选 exp，只保留映射成功的行，并把行名换成 symbol
keep      <- ens_stripped %in% names(emap)
exp2  <- exp[keep, , drop=FALSE]
rownames(exp2) <- emap[ ens_stripped[keep] ]

# 按行（基因）做 median‑centered zscore 标准化 ----
#    Z_ij = ( log2(TPM_ij+1) - median(log2(TPM_i+1)) ) / sd(log2(TPM_i+1))
log2_mat <- log2(exp2 + 1)
z_mat <- t(apply(log2_mat, 1, function(g) {
  (g - median(g)) / sd(g)
}))

xgboost$Cor <- mapply(function(lig, rec) {
  if (lig %in% rownames(z_mat) && rec %in% rownames(z_mat)) {
    cor(z_mat[lig, ], z_mat[rec, ], use="pairwise.complete.obs")
  } else {
    0
  }
}, xgboost$source, xgboost$target)
rm(exp,exp2,log2_mat,mapping,mart,my_ds)

fwrite(xgboost,"D:\\ATAC\\HCLRs\\xgboost_species_TF_cor_HCLR.txt",quote = F,sep = "\t")


##### phastCons100_30 #####
library(data.table)
xgboost <- fread("D:\\ATAC\\HCLRs\\xgboost_species_TF_cor_phastCons100_HCLR.txt",header = T,sep = "\t")
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(phastCons30way.UCSC.hg38)#phastCons100way.UCSC.hg38
library(biomaRt)

# biomaRt_Ensembl
hgnc <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# unique_genes
unique_genes <- unique(c(xgboost$source, xgboost$target))

# hgnc_symbol_entrezgene_id
gene_mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                      filters = "hgnc_symbol", 
                      values = unique_genes, 
                      mart = hgnc)

# entrezgene_id
gene_mapping$entrezgene_id <- as.character(gene_mapping$entrezgene_id)

# TxDb_phastCons
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
phast <- getGScores("phastCons30way.UCSC.hg38")#https://support.bioconductor.org/p/129276/
# phast <- phastCons100way.UCSC.hg38 #phastCons100way.UCSC.hg38
# compute_phastCons
compute_phastCons <- function(eid) {
  if (is.na(eid)) return(NA_real_)
  gr <- genes_gr[genes_gr$gene_id == eid]
  if (length(gr) == 0) return(NA_real_)
  
  # tryCatch_gscores
  scores <- tryCatch({gscores(phast, gr)}, error = function(e) {return(NULL)})
  if (is.null(scores) || length(scores) == 0) return(NA_real_)
  return(mean(scores$default, na.rm = TRUE))
}

# phastCons_vector
gene_mapping$phastCons <- sapply(gene_mapping$entrezgene_id, compute_phastCons)
phastCons_vector <- setNames(gene_mapping$phastCons, gene_mapping$hgnc_symbol)

# phastCons_vector_xgboost
xgboost$source_phastCons30 <- phastCons_vector[xgboost$source]
xgboost$target_phastCons30 <- phastCons_vector[xgboost$target]

fwrite(xgboost,"D:\\ATAC\\HCLRs\\xgboost_species_TF_cor_phastCons100_30_HCLR.txt",quote = F,sep = "\t")

##### LECIF #####
aaa <- fread("D:\\ATAC\\HCLRs\\xgboost_species_TF_cor_phastCons100_30_LECIF_HCLR.txt",sep = "\t",header = T)

##### CellChat_2species_conservatism #####
load("D:/ATAC/HCLRs/Species/CellChatDB.human.rda")
load("D:/ATAC/HCLRs/Species/CellChatDB.mouse.rda")
human <- CellChatDB.human$interaction
mouse <- CellChatDB.mouse$interaction
inter_species <- human[which(human$interaction_name%in%intersect(human$interaction_name,mouse$interaction_name)),]
rm(CellChatDB.human,CellChatDB.mouse,human,mouse)

inter_species$interaction_name <- paste0(inter_species$ligand,"--",inter_species$receptor)
inter_species <- inter_species$interaction_name

aaa$CellChat_2species_conservatism <- 0
aaa[which(aaa$pair %in% inter_species),"CellChat_2species_conservatism"] <- 1
fwrite(aaa,"D:\\ATAC\\HCLRs\\xgboost_2+3species_TF_cor_phastCons100_30_LECIF_HCLR.txt",quote = F,sep = "\t")



####### related_features_calculation #######
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(tidyr)

feat <- fread("D:\\ATAC\\HCLRs\\xgboost-model\\Features-revise\\xgboost_2+3species_TF_cor_phastCons100_30_LECIF_HCLR.txt")[,c(1:8,23,13,16,21,22)]
feat$target_LECIF[is.na(feat$target_LECIF)] <- 0
feat$source_LECIF[is.na(feat$source_LECIF)] <- 0
feat$LR_LECIF_min <- pmin(feat$target_LECIF,feat$source_LECIF) 
feat$LR_LECIF_sqmean <- sqrt(feat$target_LECIF * feat$source_LECIF)
feat$LR_LECIF_mean <- (feat$target_LECIF + feat$source_LECIF)/2

CITE_TICCOM <- fread("D:/ATAC/HCLRs/xgboost-model/code_used_data/dabiao_CITE_TICCOM_origin_195.txt",header = T,sep = "\t")[,c(1:3)]
colnames(CITE_TICCOM) <- colnames(feat)[1:3]


# Define variable sets
new_binary    <- c("CellChat_2species_conservatism", "KEGG_TF_class")
old_binary    <- c("Pathway_or_not", "LR_stimulate_inhibit", "Trans_and_Recep_either_no")
new_continuous <- c("Cor_GTEx", "LR_LECIF_min", "LR_LECIF_sqmean", "LR_LECIF_mean")
old_continuous <- c("n_references", "n_resources")

# 1. Binary intersections and proportions
bin_stats <- expand.grid(new = new_binary, old = old_binary, stringsAsFactors = FALSE) %>%
  mutate(
    count_same = map2_int(new, old, ~ sum(feat[[.x]] == feat[[.y]], na.rm = TRUE)),
    total     = nrow(feat),
    prop      = count_same / total
  )
print(bin_stats)

# 1b. Spearman correlations
cor_stats <- expand.grid(new = new_continuous, old = old_continuous, stringsAsFactors = FALSE) %>%
  mutate(
    cor_test = map2(new, old, ~ cor.test(feat[[.x]], feat[[.y]], method = "spearman", use = "complete.obs")),
    rho      = map_dbl(cor_test, ~ .x$estimate),
    p_value  = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  select(new, old, rho, p_value)
print(cor_stats)

# 2. Variance/SD for new continuous; 0/1 ratio for new binary
cont_variance <- tibble(
  var = new_continuous,
  variance = map_dbl(new_continuous, ~ var(feat[[.x]], na.rm = TRUE)),
  sd       = map_dbl(new_continuous, ~ sd(feat[[.x]], na.rm = TRUE))
)
print(cont_variance)

bin_ratio <- tibble(
  var = new_binary,
  count0 = map_int(new_binary, ~ sum(feat[[.x]] == 0, na.rm = TRUE)),
  count1 = map_int(new_binary, ~ sum(feat[[.x]] == 1, na.rm = TRUE)),
  prop0  = count0 / nrow(feat),
  prop1  = count1 / nrow(feat)
)
print(bin_ratio)

# 3. Compare CITE vs random sample
# Prepare labels
# Prepare data
cite_pairs = CITE_TICCOM$pair
feat_cite = feat[feat$pair %in% cite_pairs, ]
feat_non_cite = feat[!(feat$pair %in% cite_pairs), ]
set.seed(123)
random_sample = feat_non_cite[sample(nrow(feat_non_cite), 195), ]

# Box plot for continuous variables
for (col in new_continuous) {
  # Create plot data
  plot_data = data.frame(
    value = c(feat_cite[[col]], random_sample[[col]]),
    group = factor(rep(c("CITE_TICCOM", "Random"), each = 195), levels = c("CITE_TICCOM", "Random"))
  )
  
  # Create box plot and add p value
  p <- ggboxplot(
    plot_data, 
    x     = "group", 
    y     = "value",
    title = paste("Box plot for", col),
    xlab  = "Group", 
    ylab  = col,
    # Fill by group
    fill  = "group",
    # Custom palette: first for CITE_TICCOM, second for Random
    palette = c("#1f78b4", "#e31a1c")) +
    stat_compare_means(
      aes(label = paste0("p = ", ..p.format..)), 
      method = "wilcox.test",
      label.x  = 1.5, label.y = max(plot_data$value, na.rm = TRUE) * 1.05
    ) +
    theme_minimal() +theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  print(p)
}

# Stacked bar plot for binary variables
for (col in new_binary) {
  # Ensure factor levels are consistent
  all_levels <- unique(c(unique(feat_cite[[col]]), unique(random_sample[[col]])))
  feat_cite[[col]] <- factor(feat_cite[[col]], levels = all_levels)
  random_sample[[col]] <- factor(random_sample[[col]], levels = all_levels)
  
  # Create plot data
  plot_data_bin <- rbind(
    data.frame(group = "CITE_TICCOM", category = feat_cite[[col]], stringsAsFactors = FALSE),
    data.frame(group = "Random", category = random_sample[[col]], stringsAsFactors = FALSE)
  )
  
  # Calculate frequency
  count_table <- table(plot_data_bin$group, plot_data_bin$category)
  count_df <- as.data.frame(count_table)
  names(count_df) <- c("group", "category", "Freq")
  
  # Ensure group order
  count_df$group <- factor(count_df$group, levels = c("CITE_TICCOM", "Random"))
  
  # Create stacked bar plot
  p_bin <- ggbarplot(
    count_df, 
    x = "group", 
    y = "Freq", 
    fill = "category", 
    position = position_stack(), 
    palette = "jco",
    title = paste("Bar plot for", col),
    xlab = "Group", 
    ylab = "Count", 
    legend.title = col
  )
  
  # Fisher exact test
  fisher_res <- fisher.test(count_table)
  p_value <- fisher_res$p.value
  
  # Find max frequency to determine p value annotation position
  max_count <- max(count_df$Freq)
  # Add p value to plot
  p_bin <- p_bin + annotate(
    "text", 
    x = 1.5, 
    y = max_count * 1.1, 
    label = paste("P =", round(p_value, 4)), 
    hjust = 0.5
  )
  
  print(p_bin)
}










