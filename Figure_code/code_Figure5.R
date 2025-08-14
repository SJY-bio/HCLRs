# -----------------------------------------------------------------------------
# -----------------------------
# Step 0: Load Libraries
# -----------------------------
library(tidyverse)    # dplyr, tidyr, readr, ggplot2, etc.
library(reshape2)
library(ggplot2)
library(ggsignif)

# -----------------------------
# Step 1: Configuration - Set Paths
# -----------------------------

# Create directories
dir.create(output_dir_scores_scc,     showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_scores_mel,     showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_combined,       showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_figures,        showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Step 2: Load Metadata and Expression Data (Example)
# -----------------------------

exp <- read.csv("data/expression.csv", row.names = 1)
mela_meta <- read.csv("data/metadata.csv")

# -----------------------------
# Step 3: Compute Interaction Scores for SCC Dataset
# -----------------------------
filenames_scc <- list.files(input_dir_resources_scc, pattern = "*.txt", full.names = FALSE)

for (filename in filenames_scc) {
  filepath <- file.path(input_dir_resources_scc, filename)
  
  # Read interaction data
  scc_net_near1 <- read.table(filepath, header = TRUE, sep = "\t")
  scc_net_near1$LR <- paste(scc_net_near1$ligand, scc_net_near1$receptor, sep = "_")
  scc_net_near1$CC <- paste(scc_net_near1$source, scc_net_near1$target, sep = "_")
  
  # Skip if only Near_Near or Normal_Epithelial_Normal_Epithelial
  if (all(scc_net_near1$CC %in% c("Near_Near", "Normal_Epithelial_Normal_Epithelial"))) {
    next
  }
  
  # Filter out self-interactions
  scc_net_near_inter1 <- scc_net_near1 %>%
    filter(!CC %in% c("Near_Near", "Normal_Epithelial_Normal_Epithelial"))
  
  # Compute interaction score: sum(ligand_exp) * sum(receptor_exp)
  scc_net_near_inter1 <- scc_net_near_inter1 %>%
    mutate(LR_near_score = NA_real_)
  
  for (i in 1:nrow(scc_net_near_inter1)) {
    ligand <- scc_net_near_inter1$ligand[i]
    receptor <- scc_net_near_inter1$receptor[i]
    source_cell <- scc_net_near_inter1$source[i]
    target_cell <- scc_net_near_inter1$target[i]
  
    ligand_exp <- exp[rownames(exp) == ligand, mela_meta$cell.types == source_cell]
    receptor_exp <- exp[rownames(exp) == receptor, mela_meta$cell.types == target_cell]
    
    if (length(ligand_exp) == 0 || length(receptor_exp) == 0) {
      scc_net_near_inter1$LR_near_score[i] <- 0
    } else {
      scc_net_near_inter1$LR_near_score[i] <- sum(as.numeric(ligand_exp)) * sum(as.numeric(receptor_exp))
    }
  }
  
  # Add small random value to zero scores for log transformation
  zeros <- which(scc_net_near_inter1$LR_near_score == 0)
  if (length(zeros) > 0) {
    random_vec <- runif(length(zeros), min = 0, max = 0.001)
    scc_net_near_inter1$log_score[zeros] <- log(random_vec)
  }
  
  # Log-transform non-zero scores
  non_zeros <- which(scc_net_near_inter1$LR_near_score != 0)
  scc_net_near_inter1$log_score[non_zeros] <- log(scc_net_near_inter1$LR_near_score[non_zeros])
  
  # Normalize log score (min = 0)
  scc_net_near_inter1$log_score <- scc_net_near_inter1$log_score - min(scc_net_near_inter1$log_score)
  
  # Save scored data
  output_path <- file.path(output_dir_scores_scc, filename)
  write.table(scc_net_near_inter1, output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
# -----------------------------
# Step 4: Merge All Scored Results (Melanoma)
# -----------------------------
folder_path_mel <- output_dir_scores_mel
file_list_mel <- list.files(path = folder_path_mel, pattern = "*.txt", full.names = TRUE)

combined_data_mel <- file_list_mel %>%
  map_dfr(~ {
    dat <- read_delim(.x, delim = "\t", col_names = TRUE)
    filename <- tools::file_path_sans_ext(basename(.x))
    dat %>% mutate(resource = filename)
  })

# Save combined melanoma data
write.table(combined_data_mel,
            file.path(output_dir_combined, "all_22resource_result_melanoma.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# -----------------------------
# Step 5: Merge All Scored Results (SCC)
# -----------------------------
folder_path_scc <- output_dir_scores_scc
file_list_scc <- list.files(path = folder_path_scc, pattern = "*.txt", full.names = TRUE)

combined_data_scc <- file_list_scc %>%
  map_dfr(~ {
    dat <- read_delim(.x, delim = "\t", col_names = TRUE)
    filename <- tools::file_path_sans_ext(basename(.x))
    dat %>% mutate(resource = filename)
  })

# Rename HCLR to HCLRs for consistency
combined_data_scc$resource[combined_data_scc$resource == 'HCLR'] <- 'HCLRs'

# Save combined SCC data
write.table(combined_data_scc,
            file.path(output_dir_combined, "all_22resource_result_scc.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# -----------------------------
# Step 6: Compute Jaccard Similarity (Melanoma)
# -----------------------------
data_all_mel <- read.table(file.path(output_dir_combined, "all_22resource_result_melanoma.txt"),
                           header = TRUE, sep = "\t") %>%
  distinct(LR_pair, resource, .keep_all = TRUE)

resources_mel <- unique(data_all_mel$resource)
n_res <- length(resources_mel)
jaccard_matrix_mel <- matrix(NA, n_res, n_res, dimnames = list(resources_mel, resources_mel))

for (i in 1:n_res) {
  lr1 <- data_all_mel[data_all_mel$resource == resources_mel[i], "LR_pair"]
  for (j in 1:n_res) {
    lr2 <- data_all_mel[data_all_mel$resource == resources_mel[j], "LR_pair"]
    jaccard_matrix_mel[i, j] <- length(intersect(lr1, lr2)) / length(union(lr1, lr2))
  }
}

# Convert to long format for plotting
jaccard_long_mel <- melt(jaccard_matrix_mel)
jaccard_long_mel$value <- round(jaccard_long_mel$value, 2)

# Plot Jaccard heatmap
p_jaccard_mel <- ggplot(jaccard_long_mel, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), color = "#49525E") +
  scale_fill_gradientn(colors = c("#00CCCC", "white", "#FF8E33")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  labs(title = "Jaccard Index - Melanoma", fill = "Jaccard Index")

ggsave(p_jaccard_mel, 
       filename = file.path(output_dir_figures, "F5_a.pdf"), 
       width = 12, height = 8)

# -----------------------------
# Step 7: Plot Number of Significant LR Pairs (Bar Plot)
# -----------------------------
lr_count_mel <- data_all_mel %>%
  count(resource, name = "LR_number") %>%
  arrange(desc(LR_number))

p_bar_mel <- ggplot(lr_count_mel, aes(x = reorder(resource, -LR_number), y = LR_number, fill = resource)) +
  geom_col(width = 0.15) +
  geom_point(color = "#46505B", size = 3) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  labs(x = "", y = "LR number")

ggsave(p_bar_mel, 
       filename = file.path(output_dir_figures, "F5_a_top.pdf"), 
       width = 6, height = 4)

# -----------------------------
# Step 8: Repeat for SCC Dataset
# -----------------------------
data_all_scc <- read.table(file.path(output_dir_combined, "all_22resource_result_scc.txt"),
                           header = TRUE, sep = "\t") %>%
  distinct(LR_pair, resource, .keep_all = TRUE)

resources_scc <- unique(data_all_scc$resource)
jaccard_matrix_scc <- matrix(NA, length(resources_scc), length(resources_scc),
                             dimnames = list(resources_scc, resources_scc))

for (i in seq_along(resources_scc)) {
  lr1 <- data_all_scc[data_all_scc$resource == resources_scc[i], "LR_pair"]
  for (j in seq_along(resources_scc)) {
    lr2 <- data_all_scc[data_all_scc$resource == resources_scc[j], "LR_pair"]
    jaccard_matrix_scc[i, j] <- length(intersect(lr1, lr2)) / length(union(lr1, lr2))
  }
}

jaccard_long_scc <- melt(jaccard_matrix_scc)
jaccard_long_scc$value <- round(jaccard_long_scc$value, 2)

p_jaccard_scc <- ggplot(jaccard_long_scc, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), color = "#49525E") +
  scale_fill_gradientn(colors = c("#00CCCC", "white", "#FF8E33")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  labs(title = "Jaccard Index - SCC", fill = "Jaccard Index")

ggsave(p_jaccard_scc, 
       filename = file.path(output_dir_figures, "S5_a.pdf"), 
       width = 12, height = 8)

# Bar plot for SCC
lr_count_scc <- data_all_scc %>%
  count(resource, name = "LR_number") %>%
  arrange(desc(LR_number))

p_bar_scc <- ggplot(lr_count_scc, aes(x = reorder(resource, -LR_number), y = LR_number, fill = resource)) +
  geom_col(width = 0.15) +
  geom_point(color = "#46505B", size = 3) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  labs(x = "", y = "LR number")

ggsave(p_bar_scc, 
       filename = file.path(output_dir_figures, "S5_b.pdf"), 
       width = 7, height = 4)

# -----------------------------
# Step 9: Boxplot of log interaction scores across resources
# -----------------------------
# Melanoma
p_box_mel <- ggplot(data_all_mel, aes(x = reorder(resource, -log_score), y = log_score, fill = resource)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  labs(x = "", y = "log(interact score)")

ggsave(p_box_mel, 
       filename = file.path(output_dir_figures, "resource22_score_new_mel.pdf"), 
       width = 8, height = 8)

# SCC
p_box_scc <- ggplot(data_all_scc, aes(x = reorder(resource, log_score), y = log_score, fill = resource)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  labs(x = "", y = "log(interact score)")

ggsave(p_box_scc, 
       filename = file.path(output_dir_figures, "resource22_score_new_scc.pdf"), 
       width = 8, height = 5)

# -----------------------------
# Step 10: Wilcoxon Test for Score Differences
# -----------------------------

# -----------------------------
# Step 11: Extract Network Strength 
# -----------------------------
mel_net <- read.table(file.path(output_dir_scores_mel, "CellChatDB.txt"), 
                      sep = "\t", header = TRUE)
mel_net <- mel_net %>% select(inter_cell, log_score) %>%
  group_by(inter_cell) %>% summarise(total_score = sum(log_score)) %>%
  separate(inter_cell, into = c("source", "target"), sep = "--", remove = FALSE) %>%
  mutate(source = trimws(source), target = trimws(target))

write.table(mel_net, 
            file.path(output_dir_combined, "network", "strength_mel_CellChatDB.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

