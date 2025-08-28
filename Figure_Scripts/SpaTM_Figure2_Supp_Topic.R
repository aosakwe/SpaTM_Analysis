set.seed(1234) 
library(SpaTM)
library(scran)
library(scater) 
library(tidyverse)   
library(Matrix) 
library(ComplexHeatmap) 
library(circlize) 
library(spatialLIBD) 
library(patchwork) 
library(tictoc)
library(torch)


# Load Data
spe <- readRDS("../data/spe.rds")
#samples <- c("151673","151674","151676","151507")
#spe <- spe[,spe$sample_id %in% samples]
spe <- spe[,spe$spatialLIBD != 'Unknown']
spe$spatialLIBD <- as.factor(as.character(spe$spatialLIBD))

scenarios <- list(train = c("151673"),test = c('151674','151675','151676','151507','151508'  ,'151509',
                                               '151510',"151669",
                                               "151670",
                                               "151671",
                                               "151672"))

models <- list()
############### Number of Topics
K_vals <- c(8,10,15,25,50)
############ ###

##Remove Mitochondrial genes
library(EnsDb.Hsapiens.v75)
gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")
spe <- spe[!rownames(spe) %in% gns$gene_id,]
## Remove Ribosomal Proteins
cur_genes <- genes(EnsDb.Hsapiens.v75)$symbol
cur_genes <- cur_genes[grepl(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",cur_genes)]
genes_remove <- genes(EnsDb.Hsapiens.v75, filter = ~ symbol %in% cur_genes)$gene_id
spe <- spe[!rownames(spe) %in% genes_remove,]

spe$layer <- as.numeric(spe$spatialLIBD) - 1
n_col <- length(scenarios$test)*2 + 3
results <- matrix(0,nrow = length(K_vals),ncol = n_col)
rownames(results) <- paste(K_vals," Topics",sep = '')
colnames(results) <- c(paste(rep(scenarios$test,each = 2),
                             rep(c('',' - Smooth'),
                                 length(scenarios$test)),
                             sep = ''),'Train','Val','Runtime')
#spe <- SpatialTopicExperiment(spe,K = K)
max_iter <- 100

z <- 1
for(K in K_vals){#length(scenarios)){
  
  print(paste("Topics: ",K,sep = ''))
  spe_train <- spe[,spe$sample_id %in% scenarios$train]
  
  mlp_layers <- 1
  classes <- length(unique(spe_train$layer))
  mlp <- build_mlp(mlp_layers,K,64,classes)
  
  #Add train-validation split
  
  spe_val <- spe[,spe$sample_id == '151674']
  #Top HVGs ####
  spe_train <- computeGeometricFactors(spe_train)
  spe_train <- logNormCounts(spe_train)
  dec <- modelGeneVar(spe_train)
  top_hvgs <- getTopHVGs(dec,n = 1000)
  spe_train <- spe_train[top_hvgs,]
  spe_val <- spe_val[top_hvgs,]
  
  
  #Remove cells with no counts
  spe_train <- spe_train[,colSums(counts(spe_train)) > 0]
  spe_val <- spe_val[,colSums(counts(spe_val)) > 0]
  
  spe_train <- SpatialTopicExperiment(
    spe_train,
    guided = FALSE,
    labels = NULL,
    K,
    hvg = NULL)
  
  spe_val <- SpatialTopicExperiment(
    spe_val,
    guided = FALSE,
    labels = NULL,
    K,
    hvg = NULL)
  #Train
  D <- ncol(spe_train)
  tic()
  spe_train <- STM(spe_train,
                   'layer',
                   1,
                   max_iter,
                   TRUE,
                   TRUE,
                   FALSE,
                   thresh = 1e-10,
                   lr = 1e-05,
                   mlp = mlp,
                   mlp_layers = mlp_layers,
                   mlp_epoch = 5555500,
                   spe_val = spe_val)
  runtime <- toc()
  print("Training Complete: ")
  cat("Runtime: ")
  cat(runtime$toc)
  cat('\n')
  mlp <- metadata(spe_train)[['STM_MLP']]
  ####
  spe_train <- buildTheta(spe_train)
  spe_train <- buildPhi(spe_train)
  
  
  spe_train <-inferTopics(spe_train,1,100,FALSE,phi(spe_train))
  spe_train <- buildTheta(spe_train)
  
  
  train_pred <- pred_mlp(mlp,theta(spe_train))
  
  train_acc <- length(which(train_pred == as.numeric(spe_train$spatialLIBD)))*
    100/ncol(spe_train)
  print(paste("train Acc: ",train_acc, sep = ''))
  
  ###
  spe_val <- inferTopics(spe_val,1,100,FALSE,phi(spe_train))
  spe_val <- buildTheta(spe_val)
  
  
  val_pred <- pred_mlp(mlp,theta(spe_val))
  
  val_acc <- length(which(val_pred == as.numeric(spe_val$spatialLIBD)))*
    100/ncol(spe_val)
  print(paste("Val Acc: ",val_acc, sep = ''))
  best_phi <- phi(spe_train)
  train_res <- c(train_acc,val_acc,runtime$toc)
  ############################################
  
  #Store models
  models[[K]] <- list(spe_train,mlp)
  ############Predict Test###################
  test_res <- c()
  for(j in 1:length(scenarios$test)){
    
    spe_test <- spe[,spe$sample_id %in% scenarios$test[j]]
    ####Top HVGs ####
    spe_test <- spe_test[top_hvgs,]
    ######
    
    spe_test <- SpatialTopicExperiment(
      spe_test,
      guided = FALSE,
      labels = NULL,
      K,
      hvg = NULL)
    spe_test <- inferTopics(spe_test,1,100,FALSE,best_phi)
    spe_test <- buildTheta(spe_test)
    #Test Accuracy
    
    test_pred <- pred_mlp(mlp,theta(spe_test))
    
    
    ### RTM-based Smoothing
    nbr_list <- get_nbrs(spe_test,'sample_id','int_cell',NULL,dist = 3,loss_fun = 1)
    smooth_pred <- rtm_smooth(spe_test,test_pred,nbr_list = nbr_list)
    
    test_acc <- length(which(test_pred == as.numeric(spe_test$spatialLIBD)))*
      100/ncol(spe_test)
    print(paste(scenarios$test[j], " || Test Acc: ",test_acc, sep = ''))
    test_res <- c(test_res,test_acc)
    #results[i,j] <- c(test_acc)
    test_acc <- length(which(smooth_pred == as.numeric(spe_test$spatialLIBD)))*
      100/ncol(spe_test)
    print(paste(scenarios$test[j], " || Smooth Test Acc: ",test_acc, sep = ''))
    test_res <- c(test_res,test_acc)
    spe_test$Pred <- test_pred
    spe_test$smooth <- smooth_pred
    #vis_clus(spe_test,clustervar = 'Pred')
    #results[i,j+1] <- test_acc
    print('=============')
  }
  results[z,] <- c(test_res,train_res)
  z <- z + 1
  #####################################
} 
print(results)
write.csv(results,"./topic_fig2_benchmark_stm.csv")
