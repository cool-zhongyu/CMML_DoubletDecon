library(reticulate)
library(Matrix)
library(Seurat)
library(dplyr)
library(PRROC)
library(pbapply)
set.seed(2025)
# read python module


scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

##############################################################
# calculate doublet score on 16 benchmark datasets
##############################################################
# list to save doublet scores
score.list <- list()

# read real data from local
# change the location accordingly
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')

# loop over each dataset
for(loc in locs){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  # scrublet
  system.time({
    result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = doublet.rate, random_state = 10L)
    results <- result$scrub_doublets(min_counts=2, 
                                     min_cells=3, 
                                     min_gene_variability_pctl=85,n_prin_comps=30L)
  })
  score <- as.vector(results[[1]])
  score.list <- append(score.list, list(score))
}
# save results, change the location accordingly
saveRDS(score.list, './datasets/result/scrublet_benchmark_score.rds')


##########################################################################################################################
# pr, recall, and tnr under 10%, 20%, and 40% identification rates
##########################################################################################################################
# read score list
score.list <- readRDS('./datasets/result/scrublet_benchmark_score.rds')
# 16 benchmark datasets locations
locs <- c('./datasets/pbmc-ch.rds','./datasets/pbmc-1A-dm.rds','./datasets/pbmc-1B-dm.rds','./datasets/pbmc-1C-dm.rds')
# identification rates
rs <-c(0.1,0.2,0.4)
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
write.table(round(results,3), './datasets/result/srublet_benchmark_threshold.txt', row.names = F)

##########################################################################################################################
# pr, recall, and tnr under the thresholds determined by doubletdecon
##########################################################################################################################
# read doublet scores
score.list <- readRDS('./datasets/result/scrublet_benchmark_score.rds')
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
  data <- readRDS(loc)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  score <- score.list[[1]]
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
scrublet_res<-data.frame(ndoublet=d,precision=precisions,recall=recalls,specificity=tnrs,f1_score=f1)
write.csv(scrublet_res,"./datasets/result/scrublet_benchmark.csv",quote=F,row.names = F)
