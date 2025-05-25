library(devtools)
library(Matrix)
library(Seurat)
library(parallel)
library(DoubletDecon)
library(dplyr)
library(pbapply)

# read real data from local
# change the location accordingly
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')
res_df<-data.frame()
# loop over each dataset
for(loc in locs){
  # read data
  data <- readRDS(loc)
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
                                  rhop=1.1, 
                                  write=FALSE,
                                  useFull=FALSE, 
                                  heatmap=FALSE,
                                  centroids=TRUE,
                                  num_doubs=100, 
                                  only50=FALSE,
                                  min_uniq=4,
                                  nCores=-1)
  })
  # save prediction
  pred <- results$DRS_doublet_table$isADoublet; table(pred)
  pred <- ifelse(pred==T,1,0); table(pred)
  # result
  index.doublet <- which(label==1)
  index.singlet <- which(label==0)
  tp <- sum(pred[which(label==1)]==1); tp
  fp <- sum(pred[which(label==0)]==1); fp
  fn <- sum(pred[which(label==1)]==0); fn
  tn <- sum(pred[which(label==0)]==0); tn
  
  precision <- tp/(tp + fp); precision
  recall <- tp/(tp + fn); recall
  tnr <- tn/(tn + fp); tnr
  f1<-2*precision*recall/(precision+recall)
  
  res<-data.frame(dataset=loc,ndoublet=sum(pred[which(pred==1)]),precision=precision,recall=recall,specificity=tnr,f1_score=f1)
  res_df<-rbind(res_df,res)

}
# save the result accordingly
write.csv(res_df,"./datasets/result/benchmark_doubletdecon.csv",quote = F,row.names = F)



