# -----------------------------
# Step 1: Load and Save Data
# -----------------------------

# Read CITEdb1 and TICCom datasets
all_CITEdb1 <- read.table('data/all_CITEdb1.txt', sep = '\t', header = T, check.names = F)
colnames(all_CITEdb1)[15] <- 'source_genesymbol'
colnames(all_CITEdb1)[16] <- 'target_genesymbol'
write.table(all_CITEdb1, 'OmniPath01/CITEdb.txt', sep = '\t', quote = F, row.names = F)

all_TICCom <- read.table('data/all_TICCom.txt', sep = '\t', header = T, check.names = F)
colnames(all_TICCom)[1] <- 'source_genesymbol'
colnames(all_TICCom)[4] <- 'target_genesymbol'
write.table(all_TICCom, 'OmniPath01/TICCom.txt', sep = '\t', quote = F, row.names = F)


# -----------------------------
# Step 2: Compute Jaccard Similarity Matrix
# -----------------------------

# List all interaction files in the directory
file_dir <- "OmniPath01"
files <- list.files(file_dir, full.names = TRUE)
resource_names <- tools::file_path_sans_ext(basename(list.files(file_dir)))

# Initialize Jaccard matrices for interactions and genes
Jaccard_matrix_pair <- matrix(NA, length(files), length(files))
Jaccard_matrix_gene <- matrix(NA, length(files), length(files))
rownames(Jaccard_matrix_pair) <- resource_names
colnames(Jaccard_matrix_pair) <- resource_names
rownames(Jaccard_matrix_gene) <- resource_names
colnames(Jaccard_matrix_gene) <- resource_names

# Compute Jaccard index for each pair of resources
for (i in seq_along(files)) {
  resource1 <- read.table(files[i], stringsAsFactors = F, sep = "\t", header = T)
  interact1 <- paste(resource1$source_genesymbol, resource1$target_genesymbol, sep = "_")
  gene1 <- unique(c(resource1$source_genesymbol, resource1$target_genesymbol))
  
  for (j in seq_along(files)) {
    resource2 <- read.table(files[j], stringsAsFactors = F, sep = "\t", header = T)
    interact2 <- paste(resource2$source_genesymbol, resource2$target_genesymbol, sep = "_")
    gene2 <- unique(c(resource2$source_genesymbol, resource2$target_genesymbol))
    
    # Jaccard Index = |A ∩ B| / |A ∪ B|
    Jaccard_matrix_pair[i, j] <- length(intersect(interact1, interact2)) / length(unique(c(interact1, interact2)))
    Jaccard_matrix_gene[i, j] <- length(intersect(gene1, gene2)) / length(unique(c(gene1, gene2)))
  }
  print(paste("Finished row", i))
}

# Save Jaccard matrices
write.table(Jaccard_matrix_pair, "Jaccard_new01/Jaccard_matrix_pair.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(Jaccard_matrix_gene, "Jaccard_new01/Jaccard_matrix_gene.txt", sep = "\t", quote = F, col.names = T, row.names = T)


# -----------------------------
# Step 3: Prepare Data for Heatmap
# -----------------------------

library(reshape2)
library(ggplot2)

# Read Jaccard matrices
Jaccard_matrix_pair <- read.table("Jaccard_new01/Jaccard_matrix_pair.txt", header = T, sep = "\t", stringsAsFactors = F)
Jaccard_matrix_gene <- read.table("Jaccard_new01/Jaccard_matrix_gene.txt", header = T, sep = "\t", stringsAsFactors = F)

# Keep only lower triangle for visualization
Jaccard_matrix_pair[upper.tri(Jaccard_matrix_pair)] <- NA
Jaccard_matrix_gene[upper.tri(Jaccard_matrix_gene)] <- NA

# Melt matrices into long format
Jaccard_pair_long <- melt(as.matrix(Jaccard_matrix_pair), na.rm = TRUE)
Jaccard_gene_long <- melt(as.matrix(Jaccard_matrix_gene), na.rm = TRUE)

Jaccard_pair_long$value <- round(Jaccard_pair_long$value, 2)
Jaccard_gene_long$value <- round(Jaccard_gene_long$value, 2)

# Combine interaction and gene similarity
Jaccard_pair_long$gene_value <- Jaccard_gene_long$value


# -----------------------------
# Step 4: Load Ordering for Consistent Plot Layout
# -----------------------------

library(readxl)
order_resource <- read_xlsx("left-dotplot/lr_unique_count_per_resource.xlsx")
order_resource <- as.data.frame(order_resource)


# -----------------------------
# Step 5: Plot Jaccard Heatmap
# -----------------------------

Jaccard_plot <- ggplot(Jaccard_pair_long) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  geom_text(aes(x = Var1, y = Var2, label = value), color = "#49525E", size = 4) +
  scale_fill_gradientn(colors = c("#00CCCC", "white", "#FF8E33")) +
  geom_tile(aes(x = Var2, y = Var1, fill = gene_value)) +
  geom_text(aes(x = Var2, y = Var1, label = gene_value), color = "#49525E", size = 4) +
  scale_y_discrete(limits = rev(order_resource$Resource)) +
  scale_x_discrete(limits = order_resource$Resource) +
  xlab("") + ylab("") + labs(fill = "Jaccard Index") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20)
  )

print(Jaccard_plot)
ggsave(Jaccard_plot, 
       filename = "Jaccard_new01/Jaccard_21resource_pair_gene2.pdf", 
       width = 13, height = 9)


# -----------------------------
# Step 6: Bar Plot - LR, Ligand, Receptor Counts
# -----------------------------

resource_counts <- read_xlsx("top_barplot/resource_LRcount.xlsx")
library(dplyr)
resource_counts <- resource_counts %>% arrange(desc(Interactions_Count))

# Reshape for plotting
library(reshape2)
data_long <- melt(resource_counts, 
                  id.vars = "Resource", 
                  measure.vars = c("Receptors_Count", "Ligands_Count", "Interactions_Count"),
                  variable.name = "type", 
                  value.name = "count")

colors <- c("#4C749F", "#E7892F", "#579A4E")

bar_plot <- ggplot(data_long, aes(x = Resource, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    legend.text = element_text(size = 15)
  ) +
  labs(y = "Count", x = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = resource_counts$Resource)

print(bar_plot)
ggsave(bar_plot, "top_barplot/resource21_LR_num01.pdf", width = 14, height = 6)


# -----------------------------
# Step 7: Dot Plot - Unique LR Count per Resource
# -----------------------------

unique_LR <- read_xlsx("left-dotplot/lr_unique_count_per_resource.xlsx")
unique_LR <- as.data.frame(unique_LR)

dot_plot <- ggplot(unique_LR, aes(x = Resource, y = Unique_LR_Count_based_on_Total)) +
  geom_point(aes(color = Unique_LR_Count_based_on_Total), shape = 18, size = 6) +
  scale_color_continuous(low = "#B19429", high = "#AC7299") +
  geom_line(group = 1, color = '#5386B3', lwd = 0.5, linetype = 2) +
  scale_x_discrete(limits = rev(unique_LR$Resource)) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  ) +
  xlab("") + ylab("Unique LR count") +
  guides(color = "none")

print(dot_plot)
ggsave(dot_plot, "left-dotplot/resource_unique_LRnum.pdf", height = 7, width = 7)


# -----------------------------
# Step 8: Stacked Bar Plot - LR Overlap Proportion
# -----------------------------

# Count how many databases each LR appears in
LR_total <- data.frame(pair = character(), type = character(), stringsAsFactors = FALSE)

for (i in seq_along(files)) {
  lr <- read.table(files[i], stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  lr$pair <- paste(lr$source_genesymbol, lr$target_genesymbol, sep = "-")
  lr$type <- resource_names[i]
  LR_total <- rbind(LR_total, lr[, c("pair", "type")])
}

LR_num <- as.data.frame(table(LR_total$pair))
colnames(LR_num) <- c("pair", "num")
LR_total <- merge(LR_total, LR_num, by = "pair", all.x = TRUE)

# Compute proportion for each resource
LR_summary <- data.frame(Var1 = numeric(), Freq = numeric(), perscent = numeric(), type = character())

type_counts <- as.data.frame(table(LR_total$type))

for (i in 1:nrow(type_counts)) {
  subset_lr <- LR_total[LR_total$type == type_counts$Var1[i], ]
  freq_table <- as.data.frame(table(subset_lr$num))
  freq_table$perscent <- freq_table$Freq / nrow(subset_lr)
  freq_table$type <- type_counts$Var1[i]
  LR_summary <- rbind(LR_summary, freq_table)
}

LR_summary$type <- factor(LR_summary$type, levels = rev(unique_LR$Resource))
LR_summary$Var1 <- as.numeric(LR_summary$Var1)
LR_summary <- LR_summary[order(LR_summary$Var1), ]

stacked_plot <- ggplot(LR_summary, aes(x = type, y = perscent, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = rev(paletteer::paletteer_c("grDevices::Blue-Yellow", 30)[2:25])) +
  labs(y = "Flow Proportion (%)", x = "", fill = "Overlap Level") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

print(stacked_plot)
ggsave(stacked_plot, "left-dotplot/perscent.pdf", height = 7, width = 5)