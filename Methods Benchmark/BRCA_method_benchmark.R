############# benchmark TOP N% #############
# load library
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggh4x)

# define gold_standard and background_set
cancer_imm_CCLR <- fread("C:\\Users\\Dell\\Desktop\\TICCom_CITE_cancer_immune_CCLR.txt", header = T, sep = "\t")
HCLR <- read.table("D:/ATAC/HCLRs/xgboost-model/4301_LR_8feature.txt", header = T, sep = "\t")
gold_standard <- cancer_imm_CCLR[which(cancer_imm_CCLR$pair %in% HCLR$pair), 1:4]
colnames(gold_standard) <- c("source", "target", "ligand", "receptor")
gold_standard[] <- lapply(gold_standard, function(col) {
  if (is.character(col)) gsub(" ", "_", col) else col
})
background_set <- gold_standard %>% distinct()
aaa <- unique(paste0(background_set$source,"--",background_set$target))

# define calculate_metrics function
calculate_metrics <- function(pred_set, gold_set, background) {
  TP <- nrow(semi_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FP <- nrow(anti_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FN <- nrow(anti_join(gold_set, pred_set, by = c("source", "target", "ligand", "receptor")))
  TN <- nrow(anti_join(background, bind_rows(gold_set, pred_set), 
                       by = c("source", "target", "ligand", "receptor")))
  
  Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  Recall = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  F1 = ifelse((Precision + Recall) > 0, 
              2 * (Precision * Recall) / (Precision + Recall), 0)
  Accuracy = (TP + TN) / (TP + TN + FP + FN)
  
  return(data.frame(Precision, Recall, F1, Accuracy))
}

# define base_dir path list
base_dirs <- c(
  "D:/ATAC/HCLRs/GSE176078/CellChat",
  "D:/ATAC/HCLRs/GSE176078/CellphoneDB",
  "D:/ATAC/HCLRs/GSE176078/Nichenet",
  "D:/ATAC/HCLRs/GSE176078/celltalker",
  "D:/ATAC/HCLRs/GSE176078/iTALK",
  "D:/ATAC/HCLRs/GSE176078/commot",
  "D:/ATAC/HCLRs/GSE176078/NICHES"
)

# initialize result data frame
all_final_results <- data.frame()

# loop to process each base_dir
for (bd in base_dirs) {
  base_dir <- bd
  algorithm <- basename(base_dir)
  
  # find files
  default_files <- list.files(base_dir, pattern = "^default_.*\\.txt$", full.names = TRUE)
  hclr_files <- list.files(base_dir, pattern = "^HCLR_.*\\.txt$", full.names = TRUE)
  sample_ids <- gsub("^default_(.*)\\.txt$", "\\1", basename(default_files))
  
  # initialize result list
  results <- list(
    default = list(),
    hclr = list()
  )
  
  # process each sample
  min_count <- 1
  for (i in seq_along(sample_ids)) {
    default_data <- read.table(default_files[i], header = TRUE, sep = "\t") %>% arrange(desc(.[[5]]))
    hclr_data <- read.table(hclr_files[i], header = TRUE, sep = "\t") %>% arrange(desc(.[[5]]))
    
    if (ncol(hclr_data) > 5) {
      default_data <- default_data[which(default_data[, 6] >= min_count), 1:6]
      hclr_data <- hclr_data[which(hclr_data[, 6] >= min_count), 1:6]
    }
    
    if(nrow(default_data)!=0){
      default_data$CC <- paste0(default_data[,1],"--",default_data[,2])
      default_data <- default_data[which(default_data$CC %in% aaa),]
      default_data <- default_data[,colnames(default_data) != "CC"]
    }
    if(nrow(hclr_data)!=0){
      hclr_data$CC <- paste0(hclr_data[,1],"--",hclr_data[,2])
      hclr_data <- hclr_data[which(hclr_data$CC %in% aaa),]
      hclr_data <- hclr_data[,colnames(hclr_data) != "CC"]
    }
    
    # calculate results with different thresholds
    thresholds <- c(0.2, 0.4, 0.6, 0.8, 1.0)
    for (thresh in thresholds) {
      n_default <- round(nrow(default_data) * thresh)
      default_sub <- default_data[1:n_default, ]
      
      n_hclr <- round(nrow(hclr_data) * thresh)
      hclr_sub <- hclr_data[1:n_hclr, ]
      
      results$default[[sample_ids[i]]][[as.character(thresh)]] <- calculate_metrics(default_sub, gold_standard, background_set)
      results$hclr[[sample_ids[i]]][[as.character(thresh)]] <- calculate_metrics(hclr_sub, gold_standard, background_set)
    }
  }
  
  # calculate summary metrics
  summary_metrics <- function(method) {
    metrics_list <- c("Precision", "Recall", "F1", "Accuracy")
    summary_df <- data.frame()
    
    for (thresh in c("0.2", "0.4", "0.6", "0.8", "1")) {
      for (metric in metrics_list) {
        values <- sapply(results[[method]], function(x) x[[thresh]][[metric]])
        mean_value <- mean(values, na.rm = TRUE)
        summary_df <- rbind(summary_df, 
                            data.frame(Method = method,
                                       Threshold = as.numeric(thresh),
                                       Metric = metric,
                                       Value = mean_value))
      }
    }
    return(summary_df)
  }
  
  summary_df_default <- summary_metrics("default")
  summary_df_default$Algorithm <- algorithm
  summary_df_hclr <- summary_metrics("hclr")
  summary_df_hclr$Algorithm <- algorithm
  
  # merge current algorithm results
  algorithm_results <- bind_rows(summary_df_default, summary_df_hclr)
  
  # append to total results
  if (nrow(all_final_results) == 0) {
    all_final_results <- algorithm_results
  } else {
    all_final_results <- bind_rows(all_final_results, algorithm_results)
  }
}

# add algorithm type and sort algorithm order and metric order
all_final_results$Group <- ifelse(all_final_results$Algorithm %in% c("CellChat", "CellphoneDB", "celltalker", "iTALK", "Nichenet"), "SC", "ST")
all_final_results$Group <- factor(all_final_results$Group, levels = c("SC", "ST"))
all_final_results$Algorithm <- factor(all_final_results$Algorithm, levels = c("CellChat", "CellphoneDB", "celltalker", "iTALK", "Nichenet", "commot", "NICHES"))
all_final_results$Metric <- factor(all_final_results$Metric, levels = c("Precision", "Recall", "F1", "Accuracy"))

# generate single combined plot
p <- ggplot(all_final_results, aes(x = Threshold, y = Value, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_nested(Metric ~ Group + Algorithm, scales = "free_y", switch = "y",independent = "y", nest_line = element_line(color = "black")) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                     labels = c("20%", "40%", "60%", "80%", "100%")) +
  scale_color_manual(values = c("default" = "#ff7700", "hclr" = "#7084e7")) +
  labs(title = "", x = "Top Interactions Threshold", y = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(size = 12),
        strip.placement = "outside",
        axis.title.y = element_blank())
p
# save to PDF
# ggsave("benchmark_TOP N%_plot.pdf", p, width = 20, height = 15, units = "in")



############# benchmark TOP N #############
# load necessary library
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggh4x)

# define gold_standard and background_set
cancer_imm_CCLR <- fread("C:\\Users\\Dell\\Desktop\\TICCom_CITE_cancer_immune_CCLR.txt", header = T, sep = "\t")
HCLR <- read.table("D:/ATAC/HCLRs/xgboost-model/4301_LR_8feature.txt", header = T, sep = "\t")
gold_standard <- cancer_imm_CCLR[which(cancer_imm_CCLR$pair %in% HCLR$pair), 1:4]
colnames(gold_standard) <- c("source", "target", "ligand", "receptor")
gold_standard[] <- lapply(gold_standard, function(col) {
  if (is.character(col)) gsub(" ", "_", col) else col
})
background_set <- gold_standard %>% distinct()
aaa <- unique(paste0(background_set$source,"--",background_set$target))

# define calculate_metrics function
calculate_metrics <- function(pred_set, gold_set, background) {
  TP <- nrow(semi_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FP <- nrow(anti_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FN <- nrow(anti_join(gold_set, pred_set, by = c("source", "target", "ligand", "receptor")))
  TN <- nrow(anti_join(background, bind_rows(gold_set, pred_set), 
                       by = c("source", "target", "ligand", "receptor")))
  
  Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  Recall = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  F1 = ifelse((Precision + Recall) > 0, 
              2 * (Precision * Recall) / (Precision + Recall), 0)
  Accuracy = (TP + TN) / (TP + TN + FP + FN)
  
  return(data.frame(Precision, Recall, F1, Accuracy))
}

# define base_dir path list
base_dirs <- c(
  "D:/ATAC/HCLRs/GSE176078/CellChat",
  "D:/ATAC/HCLRs/GSE176078/CellphoneDB",
  "D:/ATAC/HCLRs/GSE176078/Nichenet",
  "D:/ATAC/HCLRs/GSE176078/celltalker",
  "D:/ATAC/HCLRs/GSE176078/iTALK",
  "D:/ATAC/HCLRs/GSE176078/commot",
  "D:/ATAC/HCLRs/GSE176078/NICHES"
)

# initialize result data frame
all_final_results <- data.frame()

# define thresholds (defined outside the loop for later use)
thresholds <- c(100, 250, 500, 1000, 2500, 5000, 10000)

# loop to process each base_dir
for (bd in base_dirs) {
  base_dir <- bd
  algorithm <- basename(base_dir)
  
  # find files
  default_files <- list.files(base_dir, pattern = "^default_.*\\.txt$", full.names = TRUE)
  hclr_files <- list.files(base_dir, pattern = "^HCLR_.*\\.txt$", full.names = TRUE)
  sample_ids <- gsub("^default_(.*)\\.txt$", "\\1", basename(default_files))
  
  # initialize result list
  results <- list(
    default = list(),
    hclr = list()
  )
  
  # process each sample
  min_count <- 1
  for (i in seq_along(sample_ids)) {
    default_data <- read.table(default_files[i], header = TRUE, sep = "\t") %>% arrange(desc(.[[5]]))
    hclr_data <- read.table(hclr_files[i], header = TRUE, sep = "\t") %>% arrange(desc(.[[5]]))
    
    if (ncol(hclr_data) > 5) {
      default_data <- default_data[which(default_data[, 6] >= min_count), 1:6]
      hclr_data <- hclr_data[which(hclr_data[, 6] >= min_count), 1:6]
    }
    
    if(nrow(default_data)!=0){
      default_data$CC <- paste0(default_data[,1],"--",default_data[,2])
      default_data <- default_data[which(default_data$CC %in% aaa),]
      default_data <- default_data[,colnames(default_data) != "CC"]
    }
    if(nrow(hclr_data)!=0){
      hclr_data$CC <- paste0(hclr_data[,1],"--",hclr_data[,2])
      hclr_data <- hclr_data[which(hclr_data$CC %in% aaa),]
      hclr_data <- hclr_data[,colnames(hclr_data) != "CC"]
    }
    
    # calculate results with different thresholds
    for (thresh in thresholds) {
      n_default <- min(thresh, nrow(default_data))
      default_sub <- default_data[1:n_default, ]
      
      n_hclr <- min(thresh, nrow(hclr_data))
      hclr_sub <- hclr_data[1:n_hclr, ]
      
      results$default[[sample_ids[i]]][[as.character(thresh)]] <- calculate_metrics(default_sub, gold_standard, background_set)
      results$hclr[[sample_ids[i]]][[as.character(thresh)]] <- calculate_metrics(hclr_sub, gold_standard, background_set)
    }
  }
  
  # calculate summary metrics
  summary_metrics <- function(method) {
    metrics_list <- c("Precision", "Recall", "F1", "Accuracy")
    summary_df <- data.frame()
    
    for (n in as.character(thresholds)) {
      for (metric in metrics_list) {
        values <- sapply(results[[method]], function(x) x[[n]][[metric]])
        mean_value <- mean(values, na.rm = TRUE)
        summary_df <- rbind(summary_df, 
                            data.frame(Method = method,
                                       N = as.numeric(n),
                                       Metric = metric,
                                       Value = mean_value))
      }
    }
    return(summary_df)
  }
  
  summary_df_default <- summary_metrics("default")
  summary_df_default$Algorithm <- algorithm
  summary_df_hclr <- summary_metrics("hclr")
  summary_df_hclr$Algorithm <- algorithm
  
  # merge current algorithm results
  algorithm_results <- bind_rows(summary_df_default, summary_df_hclr)
  
  # append to total results
  if (nrow(all_final_results) == 0) {
    all_final_results <- algorithm_results
  } else {
    all_final_results <- bind_rows(all_final_results, algorithm_results)
  }
}

# add algorithm type and sort algorithm order and metric order
all_final_results$Group <- ifelse(all_final_results$Algorithm %in% c("CellChat", "CellphoneDB", "celltalker", "iTALK", "Nichenet"), "SC", "ST")
all_final_results$Group <- factor(all_final_results$Group, levels = c("SC", "ST"))
all_final_results$Algorithm <- factor(all_final_results$Algorithm, levels = c("CellChat", "CellphoneDB", "celltalker", "iTALK", "Nichenet", "commot", "NICHES"))
all_final_results$Metric <- factor(all_final_results$Metric, levels = c("Precision", "Recall", "F1", "Accuracy"))

# generate single combined plot
p <- ggplot(all_final_results, aes(x = N, y = Value, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_nested(Metric ~ Group + Algorithm, 
               scales = "free_y", 
               switch = "y",
               independent = "y",
               nest_line = element_line(color = "black")) +
  scale_x_continuous(
    trans = "log10",  # use log scale
    breaks = thresholds,labels = scales::comma_format()  # format large numbers
  ) +
  scale_color_manual(values = c("default" = "#ff7700", "hclr" = "#7084e7")) +
  labs(title = "", x = "Top Interactions Count (log scale)", y = "") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 12),
    strip.placement = "outside",
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10,angle = 90, hjust = 1, vjust = 0.3)  # rotate X axis labels
  )
p















########## boxplot #########
# load necessary library
library(tidyverse)
library(ggpubr)
library(ggplot2)

# gold_standard_IPF <- read.table(file = "D:/ATAC/HCLRs/IPF gold standard.txt",header = T,sep = "\t")
# gold_standard_IPF$pair <- paste0(gold_standard_IPF$ligand,"--",gold_standard_IPF$receptor)
# gold_standard_CT <- read.table("D:/ATAC/HCLRs/xgboost-model/code_used_data/dabiao_CITE_TICCOM_origin_195.txt",stringsAsFactors = F,header = T,sep = "\t")[,1]
# gold_standard <- gold_standard_IPF[which(gold_standard_IPF$pair %in% gold_standard_CT),1:4]
# rownames(gold_standard) <- NULL
# background_set <- gold_standard %>% distinct()
cancer_imm_CCLR <- fread("C:\\Users\\Dell\\Desktop\\TICCom_CITE_cancer_immune_CCLR.txt",header = T,sep = "\t")
HCLR <- read.table("D:/ATAC/HCLRs/xgboost-model/4301_LR_8feature.txt",header = T,sep = "\t")
gold_standard <- cancer_imm_CCLR[which(cancer_imm_CCLR$pair %in% HCLR$pair),1:4]
colnames(gold_standard) <- c("source","target","ligand","receptor")
gold_standard[] <- lapply(gold_standard, function(col) {
  if (is.character(col)) gsub(" ", "_", col) else col
})

library(tidyverse)
library(ggplot2)
background_set <- gold_standard %>% distinct()
aaa <- unique(paste0(background_set$source,"--",background_set$target))

# define method path list
method_paths <- list(
  CellChat = "D:/ATAC/HCLRs/GSE176078/CellChat",
  CellPhoneDB = "D:/ATAC/HCLRs/GSE176078/CellPhoneDB",
  celltalker = "D:/ATAC/HCLRs/GSE176078/celltalker",
  iTALK = "D:/ATAC/HCLRs/GSE176078/iTALK/",
  Nichenet = "D:/ATAC/HCLRs/GSE176078/Nichenet",
  commot = "D:/ATAC/HCLRs/GSE176078/commot",
  NICHES = "D:/ATAC/HCLRs/GSE176078/NICHES"
)

# define calculate metrics function
calculate_metrics <- function(pred_set, gold_set, background) {
  TP <- nrow(semi_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FP <- nrow(anti_join(pred_set, gold_set, by = c("source", "target", "ligand", "receptor")))
  FN <- nrow(anti_join(gold_set, pred_set, by = c("source", "target", "ligand", "receptor")))
  
  Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  Recall = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  F1 = ifelse((Precision + Recall) > 0, 
              2 * (Precision * Recall) / (Precision + Recall), 0)
  
  return(data.frame(Precision, Recall, F1))
}

# initialize result storage
all_results <- data.frame()

# loop to process each method
for (method_name in names(method_paths)) {
  path <- method_paths[[method_name]]
  
  # get all sample files of this method
  sample_files <- list.files(path, pattern = "^HCLR_.*\\.txt$", full.names = TRUE)
  
  # extract sample ID
  sample_ids <- gsub("^HCLR_(.*)\\.txt$", "\\1", basename(sample_files))
  
  # process each sample
  for (i in seq_along(sample_files)) {
    # read data file
    sample_data <- tryCatch({
      read.table(sample_files[i], header = TRUE, sep = "\t")[, 1:5]
    }, error = function(e) {
      message(paste("Error reading file:", sample_files[i], "-", e$message))
      NULL
    })
    if(nrow(sample_data)!=0){
      sample_data$CC <- paste0(sample_data[,1],"--",sample_data[,2])
      sample_data <- sample_data[which(sample_data$CC %in% aaa),]
      sample_data <- sample_data[,colnames(sample_data) != "CC"]
    }
    
    
    if (is.null(sample_data) || nrow(sample_data) == 0) next
    
    # calculate metrics
    metrics <- calculate_metrics(
      pred_set = sample_data[, 1:4],
      gold_set = gold_standard,
      background = background_set
    )
    
    # store results
    all_results <- rbind(all_results, data.frame(
      Method = method_name,
      Sample = sample_ids[i],
      Precision = metrics$Precision,
      Recall = metrics$Recall,
      F1 = metrics$F1
    ))
  }
}

# reshape data to long format
long_results <- all_results %>%
  pivot_longer(
    cols = c(Precision, Recall, F1),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Precision", "Recall", "F1")),
    Method = factor(Method, levels = names(method_paths))
  )

# plot boxplot
gb <- ggplot(long_results, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Method", y = "Score", 
       title = "Performance Comparison of Seven Methods Using HCLR Library") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# save plot
ggsave(filename = "D:/ATAC/HCLRs/GSE176078/BRCA_method_comparison_boxplot.pdf",gb, width = 8, height = 4)


