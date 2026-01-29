####### predicted_HCLRs ##################
library(data.table)
library(Matrix)
library(patchwork)
library(pROC)
library(ggplot2)
library(xgboost)

dabiao2 <- read.table("xgboost_2+3species_TF_cor_phastCons100_30_LECIF_HCLR.txt",sep = "\t",header = T)[,c(1:8,23,13,16,21,22,9)]
dabiao2$target_LECIF[is.na(dabiao2$target_LECIF)] <- 0
dabiao2$source_LECIF[is.na(dabiao2$source_LECIF)] <- 0
dabiao2$LR_LECIF_min <- pmin(dabiao2$target_LECIF,dabiao2$source_LECIF) 
dabiao2$LR_LECIF_sqmean <- sqrt(dabiao2$target_LECIF * dabiao2$source_LECIF)
dabiao2$LR_LECIF_mean <- (dabiao2$target_LECIF + dabiao2$source_LECIF)/2
dabiao2 <- dabiao2[,c(1:13,15:17,14)]
positive <- read.table("dabiao_CITE_TICCOM_origin_195.txt",stringsAsFactors = F,header = T,sep = "\t")[,1]
negative_final <- read.table("negative1137_GO_filter_negative492.txt",stringsAsFactors = F,sep = "\t",header = T)

negative <- negative_final$LR
dabiao2$lable <- ifelse(dabiao2$pair %in% positive,1,
                        ifelse(dabiao2$pair %in% negative,0,2))

dabiao3 <- dabiao2[-which(dabiao2$lable==2),]


########### build training set #########
# traindata22 <-  Matrix(data.matrix(dabiao3[,c(4:11,15)]),sparse = T)
traindata22 <-  Matrix(data.matrix(dabiao3[,c(4:8,10:11,15)]),sparse = T)

##将自变量和因变量拼接为list##
traindata23 <- list(data=traindata22,label=dabiao3$lable)
##模型训练
set.seed(123)
xgb2 <- xgboost(data = traindata23$data,label=traindata23$label,max_depth=6, eta=0.5, objective='binary:logistic', nround=25)


#### output model importance ####
importance_matrix2 <- xgb.importance(model = xgb2)
print(importance_matrix2)
library(Ckmeans.1d.dp)

pdf("xgb2_featuer_importance.pdf",height = 5,width = 6)
xgb.ggplot.importance(importance_matrix2)
dev.off()


########4301 predicted HCLRs AUC##############
cv.res <- xgb.cv(data = traindata23$data, nfold = 10, label = traindata23$label, nrounds = 100, verbose = FALSE,
                 objective = 'binary:logistic', eval_metric = 'auc', prediction = T)
it = which.max(cv.res$evaluation_log$test_auc_mean)
best.iter = cv.res$evaluation_log$iter[it]

library(pROC)
library(ggplot2)

xgb2_ROC <- roc(response = traindata23$label,
                predictor = cv.res$pred,
                levels=c(0, 1))



ggroc(xgb2_ROC,
      color="#CC5500",
      size=1.5,
      legacy.axes = F 
)+
  theme_bw()+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        
               colour='grey', 
               linetype = 'dotdash'
  ) +
  
  #geom_point(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]]))+ 
  #geom_text(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]],label=cutOffPointText),vjust=-1)+ 
  annotate("text",x=0.6,y=0.7,label=paste("AUC = ", round(xgb2_ROC$auc,3)))

dev.off()


################################################ remove 1~7 features #####################################################
i=1
# load necessary libraries
library(xgboost)
library(Matrix)
library(pROC)
library(ggplot2)
library(paletteer)

# dabiao3_feature <- dabiao3[, c(4:11,15)]  
dabiao3_feature <- dabiao3[,c(4:8,10:11,15)]  

n_features <- ncol(dabiao3_feature)

# outer loop: k represents the number of features removed each time
for (k in 1:7) {
  # create output directory for current k value
  dir_imp <- paste0("D:/ATAC/HCLRs/xgboost-model/feature_importance/", k, "features_removed/")
  dir_roc <- paste0("D:/ATAC/HCLRs/xgboost-model/ROC/", k, "features_removed/")
  dir.create(dir_imp, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_roc, showWarnings = FALSE, recursive = TRUE)
  
  # generate all possible feature combinations (keeping n_features-k features)
  combinations <- combn(n_features, n_features - k, simplify = FALSE)
  
  # inner loop: iterate through each feature combination
  for (i in seq_along(combinations)) {
    tryCatch({
      # build training set (keeping selected features)
      selected_cols <- combinations[[i]]
      traindata22 <- Matrix(data.matrix(dabiao3_feature[, selected_cols]), sparse = TRUE)
      
      # prepare data
      traindata23 <- list(data = traindata22, label = dabiao3$lable)
      
      # train model
      xgb3 <- xgboost(
        data = traindata23$data,
        label = traindata23$label,
        max_depth = 6,
        eta = 0.5,
        objective = 'binary:logistic',
        nround = 25
      )
      
      # feature importance plot
      importance_matrix2 <- xgb.importance(model = xgb3)
      pdf(paste0(dir_imp, "combo_", i, "_feature_importance.pdf"), height = 5, width = 6)
      print(xgb.ggplot.importance(importance_matrix2))
      dev.off()
      
      # cross-validation
      cv.res <- xgb.cv(
        data = traindata23$data,
        nfold = 10,
        label = traindata23$label,
        nrounds = 100,
        verbose = FALSE,
        objective = 'binary:logistic',
        eval_metric = 'auc',
        prediction = TRUE
      )
      
      # calculate AUC
      xgb3_ROC <- roc(
        response = traindata23$label,
        predictor = cv.res$pred,
        levels = c(0, 1)
      )
      
      # get removed feature names
      removed_features <- setdiff(colnames(dabiao3_feature), colnames(dabiao3_feature)[selected_cols])
      removed_text <- paste(removed_features, collapse = ", ")
      
      # plot ROC curve
      curv_color <- paletteer_d("ggthemes::Tableau_10")[i %% 10 + 1]
      p <- ggroc(xgb3_ROC, color = curv_color, size = 1.5, legacy.axes = FALSE) +
        theme_bw() +
        geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
                     colour = 'grey', linetype = 'dotdash') +
        labs(title = paste("Removed", k, "features")) +
        annotate("text", x = 0.6, y = 0.7, 
                 label = paste0("AUC = ", round(xgb3_ROC$auc, 5), "\nRemoved: ", removed_text))
      
      pdf(paste0(dir_roc, "combo_", i, "_AUC_", round(xgb3_ROC$auc, 3), ".pdf"), 
          width = 6, height = 5)
      print(p)
      dev.off()
      
      # print progress
      cat("Completed: Removed", k, "features, Combination", i, "of", length(combinations), 
          "AUC =", round(xgb3_ROC$auc, 4), "\n")
    }, error = function(e) {
      cat("Error in k=", k, "i=", i, ":", conditionMessage(e), "\n")
    })
  }
}

############################### show AUC differences using different features ################
# set root directory
root_dir <- "D:/ATAC/HCLRs/xgboost-model/ROC"
# create empty data frame to store results
results_df <- data.frame(Feature = character(), AUC = numeric(), stringsAsFactors = FALSE)
# iterate through each feature removal folder
for (k in 1:7) {
  dir_path <- file.path(root_dir, paste0(k, "features_removed"))
  if (!dir.exists(dir_path)) next
  pdf_files <- list.files(dir_path, pattern = "\\.pdf$", full.names = TRUE)
  features_kept <- 8 - k
  feature_label <- paste0(features_kept, "_features")
  for (pdf_file in pdf_files) {
    # extract AUC value from file name
    auc_value <- str_extract(pdf_file, "AUC_([0-9\\.]+)\\.pdf")
    if (!is.na(auc_value)) {
      # extract numeric part
      auc_numeric <- as.numeric(str_replace(auc_value, "AUC_([0-9\\.]+)\\.pdf", "\\1"))
      # add to results data frame
      results_df <- rbind(results_df, data.frame(
        Feature = feature_label,
        AUC = auc_numeric,
        stringsAsFactors = FALSE
      ))
    }
  }
}
write.table(results_df,file = "D:/ATAC/HCLRs/xgboost-model/ROC/all_feature_AUC.txt",quote = F,sep = "\t",row.names = FALSE)

all_feature_AUC <- read.table("D:/ATAC/HCLRs/xgboost-model/ROC/all_feature_AUC.txt",stringsAsFactors = F,sep = "\t",header = T)
all_feature_AUC_mean <- all_feature_AUC %>%
  group_by(Feature) %>%
  summarise(
    Mean_AUC = mean(AUC),
    Min_AUC = min(AUC),
    Count = n()
  ) %>%
  arrange(desc(Mean_AUC))

aa <- ggplot(all_feature_AUC,aes(x=Feature,y=AUC,fill=Feature))+
  geom_boxplot(alpha=0.8,color="#808180",fatten=0.8,lwd=0.3,width=0.5)+
  geom_hline(yintercept = 0.98,colour="#CC5500",linetype=2,linewidth=2)+
  scale_fill_manual(values =c("#374E55","#DF8F44","#00A1D5","#B24745","#79AF97","#6A6599","#80796B"))+###人工定义颜色，color为向量
  theme_bw()+
  theme(axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12,face="plain"),
        axis.title.y=element_text(size = 16,face="plain"),
        legend.position="none")+
  ylab("AUC")+xlab("") 
aa

ggsave(aa,file="D:/ATAC/HCLRs/xgboost-model/ROC/boxplot_1-7feature.pdf",width = 8,height = 4)



###use all 687 interactions for training###
###convert independent variables to matrix##
library(Matrix)
# traindata22 <-  Matrix(data.matrix(dabiao3[,c(4:11,15)]),sparse = T)
traindata22 <-  Matrix(data.matrix(dabiao3[,c(4:8,10:11,15)]),sparse = T)

##concatenate independent variables and dependent variables as list##
traindata23 <- list(data=traindata22,label=dabiao3$lable)
##模型训练
xgb2 <- xgboost(data = traindata23$data,label=traindata23$label,max_depth=6, eta=0.5, objective='binary:logistic', nround=25)
# save(xgb2,file = "G:\\dabiao\\dabiao\\validate\\xgboost_final\\xgb2_module_all_train.Rdata")
summary(xgb2)

###===============remove complexes from 8777 interactions, leaving 7572 interactions=======================####
apply_all_pair34 <- dabiao2
###remove interactions containing complexes
apply_all_pair34_no_complex <- apply_all_pair34[which(!grepl("_",apply_all_pair34$pair)),]
apply_all_pair34_no_complex <- apply_all_pair34_no_complex[!(apply_all_pair34_no_complex$pair %in% positive),] #记得放回来

###convert independent variables to matrix##
library(Matrix)
# apply_all2 <- data.matrix(apply_all_pair34_no_complex[,c(4:11,15)])
apply_all2 <- data.matrix(apply_all_pair34_no_complex[,c(4:8,10:11,15)])

apply_all12 <- Matrix(apply_all2,sparse = T)

##concatenate independent variables and dependent variables as list##
apply_all23 <- list(data=apply_all12,label=apply_all_pair34_no_complex$location_T_P_S)

###build sparse matrix
apply_all_final2 <- xgb.DMatrix(data=apply_all23$data,label=apply_all23$label)

pre_xgb_7572 = round(predict(xgb2,newdata = apply_all_final2))

confidence_LR_index <- which(pre_xgb_7572==1)
### add CITE_TICCOM gold standard
confidence_LR <- apply_all_pair34_no_complex[confidence_LR_index,"pair"]
confidence_LR <- c(confidence_LR,positive)
# HCLRs <- read.table("D:/ATAC/HCLRs/xgboost-model/code_used_data/high_confidence_xgboost_4324.txt",stringsAsFactors = F,sep = "\t",header = T)
# length(intersect(confidence_LR,HCLRs$pair))
# length(confidence_LR)
# "GDF3--BMPR1A" %in% confidence_LR
# "CDH1--PTPRF" %in% confidence_LR
write.table(confidence_LR,file = "D:\\ATAC\\HCLRs\\xgboost-model\\4301LR.txt")
###### new feature significance ######
# posit_nega_dabiao_noloca <- dabiao3[,c(4:11,15,18)]
posit_nega_dabiao_noloca <- dabiao3[,c(4:8,10:11,15,18)]
logistic<-glm(lable~.,family="binomial",data = posit_nega_dabiao_noloca)
summary(logistic)


####### calculate F1 score for each resource #############
table(dabiao3$lable)
file <- "D:/ATAC/HCLRs/OmniPath"
dir <- dir(file)
resource <- file.path(file,dir)
resource21<-as.vector(t(sapply(dir,function(x){unlist(strsplit(x,split="[.]"))}))[,1])

F1 <- matrix(NA,21,3)
rownames(F1) <- resource21
colnames(F1) <- c("Presion","Recall","F1 score")

i=1
library(Matrix)

for(i in 1:21){
  resource_dabiao <- read.table(file=resource[i],stringsAsFactors = F,sep = "\t",header = T)
  resource_dabiao$pair <- paste(resource_dabiao$source_genesymbol,resource_dabiao$target_genesymbol,sep = "--")
  
  dabiao2_resource0 <- dabiao2[which(dabiao2$pair%in%resource_dabiao$pair),]
  # dabiao2_resource <- data.matrix(dabiao2_resource0[,c(4:11,15)])
  dabiao2_resource <- data.matrix(dabiao2_resource0[,c(4:8,10:11,15)])
  dabiao2_resource1 <- Matrix(dabiao2_resource,sparse = T)
  
  ##concatenate independent variables and dependent variables as list##
  dabiao2_resource2 <- list(data=dabiao2_resource1,label=dabiao2_resource0$location_T_P_S)
  
  ###build sparse matrix
  apply_all_final2 <- xgb.DMatrix(data=dabiao2_resource2$data,label=dabiao2_resource2$label)
  
  pre_xgb_resource21 = round(predict(xgb2,newdata = apply_all_final2))
  
  dabiao2_resource0$predict_lable <- pre_xgb_resource21
  
  F1[i,1] <- (nrow(dabiao2_resource0[which(dabiao2_resource0$lable==1 &dabiao2_resource0$predict_lable==1),])+
                nrow(dabiao2_resource0[which(dabiao2_resource0$lable==2 &dabiao2_resource0$predict_lable==1),]))/length(which(pre_xgb_resource21==1))
  F1[i,2] <- (nrow(dabiao2_resource0[which(dabiao2_resource0$lable==1 &dabiao2_resource0$predict_lable==1),])+
                nrow(dabiao2_resource0[which(dabiao2_resource0$lable==2 &dabiao2_resource0$predict_lable==1),]))/length(which(dabiao2_resource0$lable==1|dabiao2_resource0$lable==2))
  F1[i,3] <- 2* F1[i,1]*F1[i,2]/(F1[i,1]+F1[i,2])
  
  print(i)
  
}
###### predict HCLRs F1 score #####
# HCLRs <- read.table("D:/ATAC/HCLRs/xgboost-model/code_used_data/high_confidence_xgboost_4324.txt",stringsAsFactors = F,sep = "\t",header = T)

dabiao2_resource0 <- dabiao2[which(dabiao2$pair %in% confidence_LR),]
# dabiao2_resource <- data.matrix(dabiao2_resource0[,c(4:11,15)])
dabiao2_resource <- data.matrix(dabiao2_resource0[,[,c(4:8,10:11,15)])

dabiao2_resource1 <- Matrix(dabiao2_resource,sparse = T)

##concatenate independent variables and dependent variables as list##
dabiao2_resource2 <- list(data=dabiao2_resource1,lable=dabiao2_resource0$location_T_P_S)

###build sparse matrix
apply_all_final2 <- xgb.DMatrix(data=dabiao2_resource2$data,label=dabiao2_resource2$lable)

pre_xgb_resource21 = round(predict(xgb2,newdata = apply_all_final2))

dabiao2_resource0$predict_lable <- pre_xgb_resource21



precision <- (nrow(dabiao2_resource0[which(dabiao2_resource0$lable==1 &dabiao2_resource0$predict_lable==1),])+
                nrow(dabiao2_resource0[which(dabiao2_resource0$lable==2 &dabiao2_resource0$predict_lable==1),]))/length(which(pre_xgb_resource21==1))
recall <- (nrow(dabiao2_resource0[which(dabiao2_resource0$lable==1 &dabiao2_resource0$predict_lable==1),])+
             nrow(dabiao2_resource0[which(dabiao2_resource0$lable==2 &dabiao2_resource0$predict_lable==1),]))/length(which(dabiao2_resource0$lable==1|dabiao2_resource0$lable==2))
F1_1 <- 2* precision*recall/(precision+recall)

f1 <- as.matrix(data.frame(precision=precision,recall=recall,F1=F1_1))

F1 <- rbind(f1,F1)
rownames(F1)[1] <- "HCLRs"
F2 <- as.data.frame(F1)
F2$resource <- rownames(F2)

F2$scale_F1 <- scale(F2$F1)

F2$recall[2] <- 0.3990326
F2$F1[2] <- 2*F2$precision[2]*F2$recall[2]/(F2$precision[2]+F2$recall[2])
F2$scale_F1 <- scale(F2$F1)

###### plot F1 score scatter plot #######
library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
F1_22resource <- ggplot(F2,aes(resource,F1))+
  #geom_col(fill = "#BAB0AC", width = 0.1) + 
  geom_path(group=1, color = 'gray', lwd = 0.3)+
  geom_point(aes(color=F1),size=6,shape=18,stroke = 0.25)+
  #scale_colour_manual(values =myPalette(20))+
  scale_color_gradientn(colours = myPalette(20), limits=c(0.5,1))+
  theme_bw()+
  ylab("F1 score")+xlab("")+
  theme(axis.text.x = element_text(size = 12,angle =30,hjust=0.75,vjust=0.85),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15))+
  scale_x_discrete(limits=as.character(F2$resource))


















############################# HCLRs proportion in each resource  ##################################
################# see the proportion of each resource containing 4301 interactions, radar plot ###############################################################
confidence_LR <- read.table("D:\\ATAC\\HCLRs\\xgboost-model\\4301LR.txt")$x

resource22_percent <- matrix(NA,22,4)

b<-"D:/ATAC/HCLRs/OmniPath"
resource<-dir(b)
filepath<-file.path(b,resource)
resource_name <- unlist(strsplit(resource,split = ".txt"))
resource22_percent[,1] <- c(resource_name,"Union")
data21_LR_union_score <- read.table("D:\\ATAC\\HCLRs\\repair\\score_dabiao_10_col.txt",stringsAsFactors = F,sep = "\t",header = T)
#data21_LR_union_score <- dabiao_100percent[which(dabiao_100percent$pair%in%database),]##2412对，总共的8777是所有intercell互作

# HCLRs <- read.table("G:\\dabiao\\dabiao\\validate\\xgboost_final\\final_no_complex\\high_confidence_2type\\high_confidence_xgboost_4324.txt",stringsAsFactors = F,sep = "\t",header = T)


########## denominator is the number of LR in each resource ###############
i=1
for(i in 1:21){
  database <- read.table(filepath[i],stringsAsFactors = F,sep = "\t",header = T)
  data_LR <- paste(database$source_genesymbol,database$target_genesymbol,sep = "--")
  resource22_percent[i,2]<- length(which(data_LR%in%confidence_LR))/length(data_LR)
  resource22_percent[i,3] <- length(which(data_LR%in%confidence_LR))
  resource22_percent[i,4] <- length(which(data_LR%in%confidence_LR))/length(confidence_LR)
}
HCLRs_21sources_LRnum_2prop <- resource22_percent
colnames(HCLRs_21sources_LRnum_2prop) <- c("resource","proportion_in_each","LRnum","proportion_inHCLRs")
fwrite(HCLRs_21sources_LRnum_2prop,"D:\\ATAC\\HCLRs\\repair\\HCLRs_21sources_LRnum_2prop.txt",quote = F,sep = "\t")
resource22_percent <- resource22_percent[,1:3]

### add the proportion of union ##
resource22_percent[22,2]<- length(which(data21_LR_union_score$pair%in%confidence_LR))/length(data21_LR_union_score$pair)
resource22_percent[22,3] <- length(which(data21_LR_union_score$pair%in%confidence_LR))

resource22_percent[,2] <-  round(as.numeric(resource22_percent[,2]),digits=2)
resource22_percent <- as.data.frame(resource22_percent)


colnames(resource22_percent) <- c("resource","proportion","number")

############################# HCLRs proportion in each resource  ##################################
resource22_percent22 <- matrix(NA,22,3)
resource22_percent22[,1] <- c(resource_name,"Union")

i=1
for(i in 1:21){
  database <- read.table(filepath[i],stringsAsFactors = F,sep = "\t",header = T)
  data_LR <- paste(database$source_genesymbol,database$target_genesymbol,sep = "--")
  resource22_percent22[i,2]<- length(which(data_LR%in%confidence_LR))/length(confidence_LR)
  resource22_percent22[i,3] <- length(which(data_LR%in%confidence_LR))
}
### add the proportion of union ##
resource22_percent22[22,2]<- length(which(data21_LR_union_score$pair%in%confidence_LR))/length(confidence_LR)
resource22_percent22[22,3] <- length(which(data21_LR_union_score$pair%in%confidence_LR))

resource22_percent22 <- as.data.frame(resource22_percent22)
str(resource22_percent22)
resource22_percent22[,2] <-  round(as.numeric(as.character(resource22_percent22[,2])),digits=2)
resource22_percent22[,3] <-  round(as.numeric(as.character(resource22_percent22[,3])),digits=2)

colnames(resource22_percent22) <- c("resource","proportion","number")


###read resource colors####
color22 <-read.table("D:\\ATAC\\HCLRs\\repair\\21resource_color.txt",stringsAsFactors = F,sep = "\t",header = T,comment.char = "/") 
color21 <- color22[-1,]

########################### plot radar plot #####################
library(forcats) 
library(ggplot2)
resource22_percent22$proportion <- as.numeric(as.character(resource22_percent22$proportion))
resource22_percent22$number <- as.numeric(as.character(resource22_percent22$number))

resource22_percent22 <- resource22_percent22[-22,]

aa <- ggplot(resource22_percent22,aes(resource,proportion)) +                       #选择绘图列并按照分类排序
  geom_bar(aes(y=0.8,fill=resource),stat="identity",width=1, colour="white",                #先绘制整个披萨图
           alpha=0.2) +                                                                          #改变透明度
  geom_bar(stat="identity",width=1,aes(y=proportion,fill=resource),alpha=0.7) +                     #插入数值
  coord_polar() +                                                                       #变成圆形
  geom_label(aes(label=proportion,fill=resource),size=4,color="white",show.legend = FALSE)+     #为数值加标签，把'label=Per.90' 变成'label=Percentile' 以展示百分比 
  scale_fill_manual(values =color21$color2)+                                #选择披萨图的填充颜色+                                                             
  scale_y_continuous(limits = c(-0.1,1))+                                              #在中间添加白色圆圈  
  labs(fill="Resource",                                                                         #去掉图例                                        #标题
       title=NULL)+                                               #标题为资源 
  theme_minimal() +                                                                     #主题进行调整 
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust=0.5),
    plot.caption = element_text(hjust=0.5,size=6),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x=element_text(size=14,face="plain")) #设置x轴刻度标签的字体属性

table(resource21_percent_detail11$resource)

fct_reorder(resource,proportion)

ggsave(aa,file="D:\\ATAC\\HCLRs\\repair\\resource22_proportion_inHCLRs.pdf",width = 12,height = 8)


########### GO enrichment analysis for HCLRs ##############
# BiocManager::install("GOplot")
library(GOplot)
library(org.Hs.eg.db)         
library(clusterProfiler)
library(dplyr)
###read HCLRs###
# HCLRs_gene <- read.table("xgboost_final\\final_no_complex\\high_confidence_2type\\HCLRs_gene_4324_1862gene.txt",stringsAsFactors = F,sep = "\t",header = F)
# all_2type_GO <- read.table("G:\\dabiao\\dabiao\\aaaresult_figures\\Result1\\FFFinal\\GO_CC_inter_intra\\GOterm_inter_intra.txt",stringsAsFactors = F,sep = "\t",header = T)


# all_genes <- c(HCLRs_gene$V1)
all_genes <- unique(c(apply_all_pair34[which(apply_all_pair34$pair %in% confidence_LR),"target"],apply_all_pair34[which(apply_all_pair34$pair %in% confidence_LR),"source"]))
ID_genes <- bitr(all_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb =org.Hs.eg.db )

GO_BP_CC_MF <- enrichGO(gene = ID_genes$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        readable = T)
GO_HCLRs_gene <- as.data.frame(GO_BP_CC_MF)
GO_HCLRs_gene_revise <- GO_HCLRs_gene[,c(1,2,3,10,8,13)]
colnames(GO_HCLRs_gene_revise) <- c("category", "ID", "Term", "adj_pval", "zscore", "count")

###### select top 10 for BP CC MF to plot####

GO_HCLRs_gene_revise_BP <- GO_HCLRs_gene_revise[which(GO_HCLRs_gene_revise$category=="BP"),]
GO_HCLRs_gene_revise_CC <- GO_HCLRs_gene_revise[which(GO_HCLRs_gene_revise$category=="CC"),]
GO_HCLRs_gene_revise_MF <- GO_HCLRs_gene_revise[which(GO_HCLRs_gene_revise$category=="MF"),]

GO_HCLRs_gene_revise_BP10 <- GO_HCLRs_gene_revise_BP[order(GO_HCLRs_gene_revise_BP$adj_pval,decreasing = F),][1:10,]
GO_HCLRs_gene_revise_CC10 <- GO_HCLRs_gene_revise_CC[order(GO_HCLRs_gene_revise_CC$adj_pval,decreasing = F),][1:10,]
GO_HCLRs_gene_revise_MF10 <- GO_HCLRs_gene_revise_MF[order(GO_HCLRs_gene_revise_MF$adj_pval,decreasing = F),][1:10,]

GO_HCLRs_gene_revise_3term <- rbind(rbind(GO_HCLRs_gene_revise_BP10,GO_HCLRs_gene_revise_CC10),GO_HCLRs_gene_revise_MF10)

write.table(GO_HCLRs_gene_revise_3term,"D:\\ATAC\\HCLRs\\repair\\GO_HCLRs_gene_revise_top10_3term.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(GO_HCLRs_gene,"D:\\ATAC\\HCLRs\\repair\\GO_HCLRs_gene_all_3term.txt",sep = "\t",col.names = T,row.names = F,quote = F)

pdf("D:\\ATAC\\HCLRs\\repair\\GO_HCLRs_adjP_top10_3term.pdf",width = 18,height = 10)
  GOBubble(GO_HCLRs_gene_revise_3term, labels = 3,colour = c('#8064A2', '#7BA79D', '#F79646'),ID=F) #气泡图
dev.off()

#### Result2: visualize the number of intersection LR  ###
library(ggplot2)
# resource22_percent <- read.table("G:\\dabiao\\dabiao\\aaaresult_figures\\Result2\\aaresource22_union_fenbu_proportion\\HCLRs_21sources_LRnum_2prop.txt",stringsAsFactors = F,sep = "\t",header = T)
resource22_percent <- read.table("D:\\ATAC\\HCLRs\\repair\\HCLRs_21sources_LRnum_2prop.txt",stringsAsFactors = F,sep = "\t",header = T)
colnames(resource22_percent) <- c("resource","proportion_in_each","LRnum","proportion_inHCLRs")

library(paletteer)
color <- as.vector(rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30)   ))

lrnum <- ggplot(resource22_percent,aes(x=resource,y=LRnum))+
  geom_line(group = 1, color = 'gray', lwd = 0.3,linetype = 2)+
  #geom_col(fill = "#F39B7F", width = 0.08)+
  geom_point(shape=18,size=7,aes(color=proportion_in_each))+
  scale_color_gradientn(colours =color)+
  #scale_color_continuous(low = "#1F74B1", high = "#C53E01")+
  #scale_colour_manual(values =c("#E64B35","#3C5488"))+
  #scale_colour_manual(values =c("#E64B35"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 11,angle =30,hjust=0.75,vjust=0.85),
        axis.text.y = element_text(size = 11))+
  xlab("")+ylab("HCLRs count")+labs(fill="Proportion")

ggsave(lrnum,filename = "D:\\ATAC\\HCLRs\\repair\\HCLRs_21sources_LRnum_prop.pdf",width =8,height = 3)

write.table(resource22_percent,"G:\\dabiao\\dabiao\\aaaresult_figures\\Result2\\aaresource22_union_fenbu_proportion\\HCLRs_21sources_LRnum_2prop.txt",sep = "\t",col.names = T,row.names = F,quote = F)

