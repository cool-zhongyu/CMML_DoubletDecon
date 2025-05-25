library(Matrix)
library(Seurat)
library(DoubletFinder)
library(PRROC)
library(pbapply)
set.seed(2025)


# read real data from local
# change the location accordingly
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')
# loop over each dataset
system.time({
  for(loc in locs){
    i<-match(loc,locs)
    print(loc)
    data <- readRDS(loc)
    count <- data[[1]]; dim(count)
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    doublet.rate <- sum(label==1) / length(label); doublet.rate
    
    # doubletfinder
    system.time({
      ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
      doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet", min.cells = 1); doublet.seurat
      doublet.seurat <- NormalizeData(doublet.seurat)
      doublet.seurat <- ScaleData(doublet.seurat)
      doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
      doublet.seurat <- RunPCA(doublet.seurat)
      
      ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
      sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE)
      sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
      bcmvn.doublet <- find.pK(sweep.stats.doublet)
      pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
      doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = sum(label==1))
      attribute <- paste('pANN', 0.25, pK, sum(label==1), sep = '_'); attribute
    })
    score <- doublet.seurat@meta.data[[attribute]]
    score.list <- append(score.list, list(score))
    
  }
})
# save results, change the location accordingly
saveRDS(score.list, './datasets/result/doubletfinder_benchmark_score.rds')

##########################################################################################################################
# pr, recall, and tnr under 10%, 20%, and 40% identification rates
##########################################################################################################################
# read score list
score.list <- readRDS('./datasets/result/doubletfinder_benchmark_score.rds')
# 16 benchmark datasets locations
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')
# identification rates
rs <- c(0.1, 0.2, 0.4)
# result matrix; 16 rows: each dataset per row; 9 cols: precision, recall, tnr per method
results <- matrix(, nrow = length(locs), ncol = 0)

# loop over identification rates
for(r in rs){
  print('====================')
  print(r)
  precisions <- c()
  recalls <- c()
  tnrs <- c()
  f1s<-c()
  result <- matrix(data = 0, nrow = length(locs), ncol=3)
  for(i in 1:length(locs)){
    print(locs[i])
    data <- readRDS(locs[i])
    # obtain the doublet labels
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    # calculate threshold based on identification rate
    score <- score.list[[i]]
    d <- floor(length(label) * r); d
    thresh <- sort(score, decreasing = T)[d]; thresh
    # predict doublet based on threshold
    pred <- score > thresh; table(pred)
    # result
    tp <- sum(pred[which(label==1)]==1); tp
    fp <- sum(pred[which(label==0)]==1); fp
    fn <- sum(pred[which(label==1)]==0); fn
    tn <- sum(pred[which(label==0)]==0); tn
    
    precision <- tp/(tp + fp); precision
    recall <- tp/(tp + fn); recall
    tnr <- tn/(tn + fp); tnr
    f1<-2*precision*recall/(precision+recall)
    
    precisions[i] <- precision
    recalls[i] <- recall
    tnrs[i] <- tnr
    f1s[i] <- f1
  }
  result <- cbind(precisions, recalls, tnrs,f1s)
  colnames(result) <- paste(colnames(result), r, sep = '_')
  results <- cbind(results, result)
}
# changel the location and name accordingly
write.table(round(results,3), './datasets/result/doubletfinder_benchmark_threshold.txt', row.names = F)

##########################################################################################################################
# pr, recall, and tnr under the thresholds determined by doubletdecon
##########################################################################################################################
# read doublet scores
score.list <- readRDS('./datasets/result/doubletfinder_benchmark_score.rds')
# 16 benchmark datasets locations
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')
# doublet # selected by doubletdecon
d <- c(8605,1174,1361,1945)
precisions <- c()
recalls <- c()
tnrs <- c()
f1s<-c()

# loop over 16 datasets
for(i in 1:length(locs)){
  # obtain doublet labels
  data <- readRDS(locs[i])
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  score <- score.list[[i]]
  # calculate threshold based on doublet number 
  thresh <- sort(score, decreasing = T)[d[i]]
  # predict doublets
  pred <- score > thresh; table(pred)
  # result
  tp <- sum(pred[which(label==1)]==1); tp
  fp <- sum(pred[which(label==0)]==1); fp
  fn <- sum(pred[which(label==1)]==0); fn
  tn <- sum(pred[which(label==0)]==0); tn
  
  precision <- tp/(tp + fp); precision
  recall <- tp/(tp + fn); recall
  tnr <- tn/(tn + fp); tnr
  
  precisions[i] <- precision
  recalls[i] <- recall
  tnrs[i] <- tnr
}
# save the result accordingly
names(precisions) <- locs; precisions
names(recalls) <- locs; recalls
names(tnrs) <- locs; tnrs
f1<-2*precisions*recalls/(precisions+recalls);f1
doubletfinder_bench<-data.frame(ndoublet=d,precision=precisions,recall=recalls,specificity=tnrs,f1_score=f1)
write.csv(doubletfinder_bench,"./datasets/result/doubletfinder_benchmark.csv",quote=F,row.names = F)
