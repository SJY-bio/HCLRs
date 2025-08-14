
# -----------------------------
# Step 1: Load Libraries
# -----------------------------
library(ggplot2)
library(SpaTalk)
library(readxl)  # For reading .xlsx files

# -----------------------------
# Step 2: Set Paths (Use Relative Paths)
# -----------------------------
data_dir <- "data/raw/"
results_dir_near <- "results/interaction/Near_fov_union/"
results_dir_far  <- "results/interaction/Far_fov_union/"
figures_dir      <- "results/figures/"

# Create directories
dir.create(results_dir_near, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir_far,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,      showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Step 3: Load Data
# -----------------------------
cell_meta_lung13 <- read.table(
  file.path(data_dir, "basic_Lung13_data.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

Lung13_xy_meta <- read.table(
  file.path(data_dir, "Lung13_xy_meta.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

Lung13_exp_final <- read.table(
  file.path(data_dir, "Lung13_exp.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# High-confidence ligand-receptor pairs
high_conf_LR <- read.table(
  file.path(data_dir, "xgboost_confidence_LR_01_8777.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

High_LR_pairs <- data.frame(
  ligand = high_conf_LR$source,
  receptor = high_conf_LR$target,
  species = "Human"
)

# Cell-cell distance info (tumor vs others)
CC_dist_tumor_other <- read.table(
  file.path(data_dir, "CC_dist_tumor_other_meta_near_far.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# -----------------------------
# Step 4: Define Helper Function to Process Each FOV
# -----------------------------
process_fov <- function(i, meta_data, exp_data, dist_data, lr_pairs, 
                        results_near_dir, results_far_dir) {
  
  message("Processing FOV ", i)
  
  # Extract cells in current FOV
  meta_20 <- meta_data[meta_data$fov == i, , drop = FALSE]
  
  # Filter CC distance data for current FOV
  dist_fov <- dist_data[
    dist_data$to %in% meta_20$cell_ID & 
    dist_data$from %in% meta_20$cell_ID, 
  ]
  
  # Sort by distance and classify into Near/Far/Other
  dist_fov <- dist_fov[order(dist_fov$distance), ]
  n_total <- nrow(dist_fov)
  n_near <- round(n_total / 4)
  n_far  <- round(n_total * 3 / 4)
  
  dist_fov$near_far_other <- "Other"
  dist_fov$near_far_other[1:n_near] <- "Near"
  dist_fov$near_far_other[(n_far + 1):n_total] <- "Far"
  
  # Extract unique cells in Near and Far groups
  extract_cell_types <- function(df, group) {
    from_cells <- df[df$near_far_other == group, c("from", "from_ct")]
    to_cells   <- df[df$near_far_other == group, c("to", "to_ct")]
    colnames(from_cells) <- colnames(to_cells) <- c("cell", "celltype")
    unique(rbind(from_cells, to_cells))
  }
  
  near_cells <- extract_cell_types(dist_fov, "Near")
  far_cells  <- extract_cell_types(dist_fov, "Far")
  
  # Add spatial coordinates
  near_with_xy <- merge(near_cells, meta_20, by.x = "cell", by.y = "cell_ID")
  far_with_xy  <- merge(far_cells,  meta_20, by.x = "cell", by.y = "cell_ID")
  
  # Extract expression data for Near/Far cells
  exp_near <- exp_data[, colnames(exp_data) %in% near_with_xy$cell, drop = FALSE]
  exp_far  <- exp_data[, colnames(exp_data) %in% far_with_xy$cell,  drop = FALSE]
  
  # Reorder metadata to match expression matrix
  near_with_xy <- near_with_xy[match(colnames(exp_near), near_with_xy$cell), ]
  far_with_xy  <- far_with_xy[match(colnames(exp_far),  far_with_xy$cell), ]
  
  # Create SpaTalk metadata
  meta_near <- data.frame(
    cell = near_with_xy$cell,
    x = near_with_xy$sdimx,
    y = near_with_xy$sdimy
  )
  
  meta_far <- data.frame(
    cell = far_with_xy$cell,
    x = far_with_xy$sdimx,
    y = far_with_xy$sdimy
  )
  
  # Build SpaTalk objects
  obj_near <- createSpaTalk(
    st_data = exp_near,
    st_meta = meta_near,
    if_st_is_sc = TRUE,
    spot_max_cell = 1,
    celltype = near_with_xy$celltype,
    species = "Human"
  )
  
  obj_far <- createSpaTalk(
    st_data = exp_far,
    st_meta = meta_far,
    if_st_is_sc = TRUE,
    spot_max_cell = 1,
    celltype = far_with_xy$celltype,
    species = "Human"
  )
  
  # Find LR pathways
  obj_near <- find_lr_path(object = obj_near, lrpairs = lr_pairs, pathways = pathways)
  obj_far  <- find_lr_path(object = obj_far,  lrpairs = lr_pairs, pathways = pathways)
  
  # Define tumor-other interaction pairs
  ct_types_near <- unique(obj_near@meta[["rawmeta"]][["celltype"]])
  tumor_other_near <- data.frame(
    sender = "tumor",
    receiver = ct_types_near[ct_types_near != "tumor"]
  )
  
  ct_types_far <- unique(obj_far@meta[["rawmeta"]][["celltype"]])
  tumor_other_far <- data.frame(
    sender = "tumor",
    receiver = ct_types_far[ct_types_far != "tumor"]
  )
  
  # Predict interactions: tumor → other and other → tumor
  for (j in 1:nrow(tumor_other_near)) {
    sender   <- tumor_other_near$sender[j]
    receiver <- tumor_other_near$receiver[j]
    
    # tumor → other
    result1 <- tryCatch({
      obj_temp <- dec_cci(obj_near, celltype_sender = sender, celltype_receiver = receiver, n_neighbor = 2)
      save(obj_temp, file = file.path(results_near_dir, paste0("fov_", i, "_tumor_", receiver, "_near1.Rdata")))
      obj_temp@lrpair
    }, error = function(e) NULL)
    
    # other → tumor
    result2 <- tryCatch({
      obj_temp <- dec_cci(obj_near, celltype_sender = receiver, celltype_receiver = sender, n_neighbor = 2)
      save(obj_temp, file = file.path(results_near_dir, paste0("fov_", i, "_tumor_", receiver, "_near2.Rdata")))
      obj_temp@lrpair
    }, error = function(e) NULL)
    
    # Combine results
    combined_result <- na.omit(bind_rows(result1, result2))
    if (nrow(combined_result) > 0) {
      combined_result$near_far <- "Near"
      combined_result$fov <- paste0("fov", i)
      combined_result$CC_inter <- paste0("tumor--", receiver)
      combined_result$resource <- "HCLRs"
      
      write.table(
        combined_result,
        file = file.path(results_near_dir, paste0("fov_", i, "_tumor_", receiver, "_near.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
      )
    }
  }
  
  for (k in 1:nrow(tumor_other_far)) {
    sender   <- tumor_other_far$sender[k]
    receiver <- tumor_other_far$receiver[k]
    
    # tumor → other
    result1 <- tryCatch({
      obj_temp <- dec_cci(obj_far, celltype_sender = sender, celltype_receiver = receiver, n_neighbor = 2)
      save(obj_temp, file = file.path(results_far_dir, paste0("fov_", i, "_tumor_", receiver, "_far1.Rdata")))
      obj_temp@lrpair
    }, error = function(e) NULL)
    
    # other → tumor
    result2 <- tryCatch({
      obj_temp <- dec_cci(obj_far, celltype_sender = receiver, celltype_receiver = sender, n_neighbor = 2)
      save(obj_temp, file = file.path(results_far_dir, paste0("fov_", i, "_tumor_", receiver, "_far2.Rdata")))
      obj_temp@lrpair
    }, error = function(e) NULL)
    
    combined_result <- na.omit(bind_rows(result1, result2))
    if (nrow(combined_result) > 0) {
      combined_result$near_far <- "Far"
      combined_result$fov <- paste0("fov", i)
      combined_result$CC_inter <- paste0("tumor--", receiver)
      combined_result$resource <- "HCLRs"
      
      write.table(
        combined_result,
        file = file.path(results_far_dir, paste0("fov_", i, "_tumor_", receiver, "_far.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
      )
    }
  }
  
  message("FOV ", i, " completed.")
}

# -----------------------------
# Step 5: Loop Over All FOVs (1 to 20)
# -----------------------------
n_fovs <- 20

for (fov_index in 1:n_fovs) {
  process_fov(
    i = fov_index,
    meta_data = Lung13_xy_meta,
    exp_data = Lung13_exp_final,
    dist_data = CC_dist_tumor_other,
    lr_pairs = High_LR_pairs,
    results_near_dir = results_dir_near,
    results_far_dir = results_dir_far
  )
}

# -----------------------------
# Step 6: Visualize Cell Type Proportion Across FOVs
# -----------------------------
# Load proportion data (assumes this file exists)
fov20_near_ctnum <- read_excel("data/raw/all22resource_20fov_LR_ct_num.xlsx")

# Load color scheme
color_data <- read.table(
  file.path(data_dir, "21resource_color.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE,
  comment.char = "$"
)

# Count unique cell types per FOV (excluding tumor)
num_file_dir <- "data/raw/fov_near_far_ct_num"
file_list <- list.files(num_file_dir, full.names = TRUE)

ct_count_per_fov <- sapply(1:n_fovs, function(i) {
  allnum <- read.table(file_list[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  n_non_tumor <- sum(!duplicated(allnum[allnum$near_far == "Near", "cell_ID"])) - 1  # minus tumor
  return(n_non_tumor)
})

fov_summary <- data.frame(fov = paste0("fov", 1:n_fovs), num = ct_count_per_fov)

# Filter Near data and merge
near_data <- fov20_near_ctnum[fov20_near_ctnum$Distance == "Near", ]
colnames(near_data)[1] <- "fov"
merged_data <- merge(fov_summary, near_data, by = "fov")
merged_data$prop <- as.numeric(merged_data$cell_type_count) / as.numeric(merged_data$num)

# Convert fov to factor for ordered plotting
merged_data$fov <- factor(merged_data$fov, levels = paste0("fov", 1:20))

# Plot proportions
p <- ggplot(merged_data, aes(x = fov, y = prop, color = resource)) +
  geom_line(aes(group = resource), linewidth = 0.6) +
  geom_point(size = 3, shape = 16) +
  scale_color_manual(values = setNames(color_data$color2, color_data$resource)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, hjust = 0.75, vjust = 0.85),
    axis.text.y = element_text(size = 12)
  ) +
  labs(x = "", y = "Proportion")

# Save plot
ggsave(p, file = file.path(figures_dir, "cell_type_proportion.pdf"), width = 11, height = 5)

# Save data
write.table(
  merged_data,
  file.path(figures_dir, "cell_type_proportion.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

message("All FOVs processed and visualization saved.")