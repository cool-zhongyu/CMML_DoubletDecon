library(devtools)
library(Matrix)
library(Seurat)
library(parallel)
library(DoubletDecon)
library(dplyr)
library(pbapply)
# read real data from local
# change the location accordingly
locs <- c('./datasets/pbmc-ch.rds')
# read data
data <- readRDS(locs)
count <- data[[1]]; dim(count)
label <- data[[2]]; table(label)
label <- ifelse(label == 'doublet', 1, 0); table(label)

# preprocess
doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet",
                                     min.cells = 1 , min.features = 1); doublet.seurat
doublet.seurat <- NormalizeData(doublet.seurat)
doublet.seurat <- ScaleData(doublet.seurat)
doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
doublet.seurat <- RunPCA(doublet.seurat)
doublet.seurat <- FindNeighbors(doublet.seurat, dims = 1:10)
doublet.seurat <- FindClusters(doublet.seurat, resolution = 0.5)
doublet.seurat <- RunUMAP(doublet.seurat,dims = 1:10)

rhops<-c(0.9,1.0,1.1)
res_df<-data.frame()

# loop over each dataset
for(rhop in rhops){

  
  # doubletdecon
  system.time({
    newFiles <- Improved_Seurat_Pre_Process(doublet.seurat, num_genes=50, write_files=FALSE)
    results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                                  groupsFile=structure(newFiles$newGroupsFile,class="matrix"), 
                                  filename='DoubletDecon', 
                                  location='./datasets/result/',
                                  fullDataFile=NULL, 
                                  removeCC=FALSE, 
                                  species="hsa", 
                                  rhop=rhop, 
                                  write=FALSE, 
                                  PMF=TRUE, 
                                  useFull=FALSE, 
                                  heatmap=FALSE,
                                  centroids=TRUE,
                                  num_doubs=100, 
                                  only50=FALSE,
                                  min_uniq=4,
                                  nCores=-1)

  })
  # save prediction
  pred_rescue <- rownames(results$DRS_doublet_table) %in% rownames(results$Final_doublets_groups)
  pred<-results$DRS_doublet_table$isADoublet
  pred_rescue <- ifelse(pred_rescue==T,1,0); table(pred_rescue)
  pred <- ifelse(pred==T,1,0); table(pred)
  # result
  index.doublet <- which(label==1)
  index.singlet <- which(label==0)
  tp <- sum(pred[which(label==1)]==1); tp
  fp <- sum(pred[which(label==0)]==1); fp
  fn <- sum(pred[which(label==1)]==0); fn
  tn <- sum(pred[which(label==0)]==0); tn
  
  tp_rescue <- sum(pred_rescue[which(label==1)]==1); tp_rescue
  fp_rescue <- sum(pred_rescue[which(label==0)]==1); fp_rescue
  fn_rescue <- sum(pred_rescue[which(label==1)]==0); fn_rescue
  tn_rescue <- sum(pred_rescue[which(label==0)]==0); tn_rescue
  
  precision <- tp/(tp + fp); precision
  recall <- tp/(tp + fn); recall
  tnr <- tn/(tn + fp); tnr
  
  precision_rescue <- tp_rescue/(tp_rescue + fp_rescue); precision_rescue
  recall_rescue <- tp_rescue/(tp_rescue + fn_rescue); recall_rescue
  tnr_rescue <- tn_rescue/(tn_rescue + fp_rescue); tnr_rescue
  
  res<-data.frame(rhop=rhop,rescue=FALSE,ndoublet=sum(pred[which(pred==1)]),precision=precision,recall=recall,specificity=tnr)
  res_rescue<-data.frame(rhop=rhop,rescue=TRUE,ndoublet=sum(pred_rescue[which(pred_rescue==1)]),precision=precision_rescue,recall=recall_rescue,specificity=tnr_rescue)
  
  res_df<-rbind(res_df,res)
  res_df<-rbind(res_df,res_rescue)

}
res_df$f1<-2*res_df$precision*res_df$recall/(res_df$precision+res_df$recall)
# save the result accordingly
write.csv(round(res_df,3),"./datasets/result/doubletdecon.csv",quote = F,row.names = F)
