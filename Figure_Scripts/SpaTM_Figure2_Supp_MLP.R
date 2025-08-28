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
library(EnsDb.Hsapiens.v75)

# Load Data
spe <- readRDS("../data/spe.rds")
#samples <- c("151673","151674","151676","151507")
#spe <- spe[, spe$sample_id %in% samples]
spe <- spe[, spe$spatialLIBD != 'Unknown']
spe$spatialLIBD <- as.factor(as.character(spe$spatialLIBD))

# Define scenarios (train, validation, test)
scenarios <- list(train = c("151673"), test = c('151674','151675','151676','151507','151508'  ,'151509',
                                                '151510',"151669",
                                                "151670",
                                                "151671",
                                                "151672"))

# Fixed number of topics and maximum iterations
K <- 8
max_iter <- 100

# Remove Mitochondrial genes
gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")
spe <- spe[!rownames(spe) %in% gns$gene_id,]

# Remove Ribosomal Proteins
cur_genes <- genes(EnsDb.Hsapiens.v75)$symbol
cur_genes <- cur_genes[grepl(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA", cur_genes)]
genes_remove <- genes(EnsDb.Hsapiens.v75, filter = ~ symbol %in% cur_genes)$gene_id
spe <- spe[!rownames(spe) %in% genes_remove,]

# Set layer (assumes layers are numeric: subtract one)
spe$layer <- as.numeric(spe$spatialLIBD) - 1

# Define the MLP structures to explore.
# Use numeric(0) for no hidden layer.
mlp_structs <- list(
  "none" = numeric(0),
  "16" = c(16),
  "32" = c(32),
  "64" = c(64),
  "64_64" = c(64,64),
  "64_64_64" = c(64,64,64)
)

# Preallocate results matrix.
# Columns: Within-Sample, Within-Sample (Smooth), Cross-Sample, Cross-Sample (Smooth), Train, Val, Runtime
n_col <- length(scenarios$test)*2 + 3
results <- matrix(0, nrow = length(mlp_structs), ncol = n_col)
rownames(results) <- names(mlp_structs)
colnames(results) <-  c(paste(rep(scenarios$test,each = 2),
                              rep(c('',' - Smooth'),
                                  length(scenarios$test)),
                              sep = ''),'Train','Val','Runtime')
model_list <- list()
row_index <- 1

for(struct_name in names(mlp_structs)) {
  cat("Evaluating MLP structure:", struct_name, "\n")
  
  ## Prepare training and validation data
  spe_train <- spe[, spe$sample_id %in% scenarios$train]
  # Use sample 151674 for validation
  spe_val <- spe[, spe$sample_id == "151674"]
  
  # Compute top HVGs (forcing exactly 1000)
  spe_train <- computeGeometricFactors(spe_train)
  spe_train <- logNormCounts(spe_train)
  dec <- modelGeneVar(spe_train)
  top_hvgs <- getTopHVGs(dec, n = 1000)
  spe_train <- spe_train[top_hvgs, ]
  spe_val <- spe_val[top_hvgs, ]
  
  # Remove cells with no counts
  spe_train <- spe_train[, colSums(counts(spe_train)) > 0]
  spe_val <- spe_val[, colSums(counts(spe_val)) > 0]
  
  # Create SpatialTopicExperiment objects for train and validation
  spe_train <- SpatialTopicExperiment(
    spe_train,
    guided = FALSE,
    labels = NULL,
    K = K,
    hvg = NULL
  )
  spe_val <- SpatialTopicExperiment(
    spe_val,
    guided = FALSE,
    labels = NULL,
    K = K,
    hvg = NULL
  )
  
  # Determine classes from training data
  classes <- length(unique(spe_train$layer))
  
  # Build the MLP using the current structure.
  current_structure <- mlp_structs[[struct_name]]
  # Here we assume build_mlp accepts a vector of hidden units.
  if(length(current_structure) == 0) {
    mlp <- build_mlp(0, K, 64, classes)  # call for no hidden layer
  } else {
    mlp <- build_mlp(length(current_structure), K, current_structure, classes)
  }
  
  # Train the STM model with the specified MLP, tracking runtime.
  tic()
  spe_train <- STM(spe_train,
                   'layer',
                   1,  # number of iterations for topic refinement
                   max_iter,
                   TRUE,  # other TRUE/FALSE flags as in your original script
                   TRUE,
                   FALSE,
                   thresh = 1e-10,
                   lr = 1e-05,
                   mlp = mlp,
                   mlp_layers = if(length(current_structure)==0) 0 else length(current_structure),
                   mlp_epoch = 100,
                   spe_val = spe_val)
  runtime <- toc(quiet = TRUE)
  cat("Training complete for structure:", struct_name, "\n")
  cat("Runtime:", runtime$toc, "\n")
  
  # Get trained MLP from metadata.
  mlp <- metadata(spe_train)[['STM_MLP']]
  
  # Build Theta and Phi then re-infer topics on training data.
  spe_train <- buildTheta(spe_train)
  spe_train <- buildPhi(spe_train)
  spe_train <- inferTopics(spe_train, 1, 100, FALSE, phi(spe_train))
  spe_train <- buildTheta(spe_train)
  
  # Prediction on training data.
  train_pred <- pred_mlp(mlp, theta(spe_train))
  train_acc <- length(which(train_pred == as.numeric(spe_train$spatialLIBD))) * 100 / ncol(spe_train)
  
  # Prediction on validation data.
  spe_val <- inferTopics(spe_val, 1, 100, FALSE, phi(spe_train))
  spe_val <- buildTheta(spe_val)
  val_pred <- pred_mlp(mlp, theta(spe_val))
  val_acc <- length(which(val_pred == as.numeric(spe_val$spatialLIBD))) * 100 / ncol(spe_val)
  
  # Save best phi from training to use for test inference.
  best_phi <- phi(spe_train)
  
  # Evaluate on test data.
  test_results <- c()
  for(test_sample in scenarios$test) {
    spe_test <- spe[, spe$sample_id %in% test_sample]
    spe_test <- spe_test[top_hvgs, ]
    spe_test <- SpatialTopicExperiment(
      spe_test,
      guided = FALSE,
      labels = NULL,
      K = K,
      hvg = NULL
    )
    spe_test <- inferTopics(spe_test, 1, 100, FALSE, best_phi)
    spe_test <- buildTheta(spe_test)
    
    # Test prediction
    test_pred <- pred_mlp(mlp, theta(spe_test))
    test_acc <- length(which(test_pred == as.numeric(spe_test$spatialLIBD))) * 100 / ncol(spe_test)
    cat(test_sample, "|| Test Acc:", test_acc, "\n")
    test_results <- c(test_results, test_acc)
    
    # RTM-based smoothing for test prediction.
    nbr_list <- get_nbrs(spe_test, 'sample_id', 'int_cell', NULL, dist = 3, loss_fun = 1)
    smooth_pred <- rtm_smooth(spe_test, test_pred, nbr_list = nbr_list)
    smooth_acc <- length(which(smooth_pred == as.numeric(spe_test$spatialLIBD))) * 100 / ncol(spe_test)
    cat(test_sample, "|| Smooth Test Acc:", smooth_acc, "\n")
    test_results <- c(test_results, smooth_acc)
    cat("=============\n")
  }
  
  # Record the metrics in the following order:
  # Within-Sample, Within-Sample (Smooth), Cross-Sample, Cross-Sample (Smooth), Train, Val, Runtime
  results[row_index, ] <- c(test_results, train_acc, val_acc, runtime$toc)
  
  model_list[[struct_name]] <- list(model = mlp, spe_train = spe_train)
  row_index <- row_index + 1
}

print(results)
write.csv(results, "./mlp_benchmark_stm.csv")