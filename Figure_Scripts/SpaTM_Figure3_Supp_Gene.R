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
library(ape)  # For Moran.I

# Load Data
spe <- readRDS("./spe.rds")
samples <- c("151673","151674","151676","151507")
spe <- spe[, spe$sample_id %in% samples]
spe <- spe[, spe$spatialLIBD != 'Unknown']
spe$spatialLIBD <- as.factor(as.character(spe$spatialLIBD))

scenarios <- list(train = c("151673"), test = c('151676','151507'))

models <- list()
############### Number of Topics and Gene Set Sizes
gene_vals <- c(500,1000,5000,10000,25000)
K <- 10

## Remove Mitochondrial genes
gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")
spe <- spe[!rownames(spe) %in% gns$gene_id,]
## Remove Ribosomal Proteins
cur_genes <- genes(EnsDb.Hsapiens.v75)$symbol
cur_genes <- cur_genes[grepl(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA", cur_genes)]
genes_remove <- genes(EnsDb.Hsapiens.v75, filter = ~ symbol %in% cur_genes)$gene_id
spe <- spe[!rownames(spe) %in% genes_remove,]

spe$layer <- as.numeric(spe$spatialLIBD) - 1

# Prepare results matrix (columns: Within-Sample & Cross-Sample average SpaTM-S Moran's I)
results <- matrix(0, nrow = length(gene_vals), ncol = 2)
rownames(results) <- paste(gene_vals, " Genes", sep = '')
colnames(results) <- c('Within-Sample','Cross-Sample')

# Lists to store per-gene Moran's I values for each gene set
all_moran676 <- list()
all_moran507 <- list()

max_iter <- 100
z <- 1
for(n_gene in gene_vals){  # Iterate over gene set sizes
  
  cat("NGenes:", n_gene, "\n")
  spe_train <- spe[, spe$sample_id %in% scenarios$train]
  
  mlp_layers <- 1
  classes <- length(unique(spe_train$layer))
  mlp <- build_mlp(mlp_layers, K, 64, classes)
  
  # Add train-validation split
  spe_val <- spe[, spe$sample_id == '151674']
  
  # Select top HVGs (exactly n_gene)
  spe_train <- computeGeometricFactors(spe_train)
  spe_train <- logNormCounts(spe_train)
  dec <- modelGeneVar(spe_train)
  top_hvgs <- getTopHVGs(dec, n = n_gene)
  cat("Top HVGs found:", length(top_hvgs), "\n")
  spe_train <- spe_train[top_hvgs, ]
  spe_val <- spe_val[top_hvgs, ]
  
  # Remove cells with no counts
  spe_train <- spe_train[, colSums(counts(spe_train)) > 0]
  spe_val <- spe_val[, colSums(counts(spe_val)) > 0]
  
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
  
  # Train STM model
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
                   mlp_epoch = 500,
                   spe_val = spe_val)
  runtime <- toc()
  cat("Training Complete; Runtime:", runtime$toc, "\n")
  mlp <- metadata(spe_train)[['STM_MLP']]
  
  spe_train <- buildTheta(spe_train)
  spe_train <- buildPhi(spe_train)
  
  best_phi <- phi(spe_train)
  
  ############################################
  # Predict on test data by combining the two test samples
  spe_test <- spe[, spe$sample_id %in% scenarios$test]
  spe_test <- spe_test[top_hvgs, ]
  spe_test <- SpatialTopicExperiment(
                spe_test,
                guided = FALSE,
                labels = NULL,
                K,
                hvg = NULL)
  spe_test <- inferTopics(spe_test, 1, 100, FALSE, best_phi)
  spe_test <- buildTheta(spe_test)
  
  # Reconstruct gene-expression from STM:
  # stm_gex = phi * t(theta) with columnwise normalization
  stm_phi <- phi(spe_train)
  stm_gex <- stm_phi %*% t(theta(spe_test))
  stm_gex <- apply(stm_gex, 2, function(a) { a/sum(a) })
  
  # Add total counts and create a new assay "stm_gex"
  spe_test$TotalCount <- colSums(counts(spe_test))
  assay(spe_test, 'stm_gex') <- sweep(stm_gex, 2, spe_test$TotalCount, `*`)
  
  # Filter DE genes using DLPFC markers (only keep genes in both spe_test and marker list)
  dlpfc_markers <- readRDS("./DLPFC_markers.rds")
  dlpfc_markers <- intersect(rownames(spe_test), dlpfc_markers$gene)
  cat("DE genes (DLPFC markers) after intersect:", length(dlpfc_markers), "\n")
  
  # Separate test data by sample
  spe676 <- spe_test[, spe_test$sample_id == '151676']
  spe507 <- spe_test[, spe_test$sample_id == '151507']
  
  # Normalize the reconstructed expression (stm_gex assay) for each sample
  spe676 <- logNormCounts(spe676,
                          name = 'stm_lognorm',
                          assay.type = 'stm_gex',
                          size.factors = NULL)
  spe507 <- logNormCounts(spe507,
                          name = 'stm_lognorm',
                          assay.type = 'stm_gex',
                          size.factors = NULL)
  
  # Compute inverse distances for Moran's I calculation (ensure array_row and array_col exist)
  dist676 <- SpaTM:::get_dist_cpp(as.matrix(colData(spe676)[, c('array_row','array_col')]))
  dist676 <- 1/dist676
  diag(dist676) <- 0
  
  dist507 <- SpaTM:::get_dist_cpp(as.matrix(colData(spe507)[, c('array_row','array_col')]))
  dist507 <- 1/dist507
  diag(dist507) <- 0
  
  # Compute Moran's I for each DE gene
  moran_res676 <- list()
  for(cur_gene in dlpfc_markers){
    m_uncorr <- Moran.I(logcounts(spe676)[cur_gene, ], dist676)$observe
    m_corr <- Moran.I(assay(spe676, "stm_lognorm")[cur_gene, ], dist676)$observed
    moran_res676[[cur_gene]] <- c(Uncorrected = m_uncorr, SpaTM_S = m_corr)
  }
  
  moran_res507 <- list()
  for(cur_gene in dlpfc_markers){
    m_uncorr <- Moran.I(logcounts(spe507)[cur_gene, ], dist507)$observe
    m_corr <- Moran.I(assay(spe507, "stm_lognorm")[cur_gene, ], dist507)$observed
    moran_res507[[cur_gene]] <- c(Uncorrected = m_uncorr, SpaTM_S = m_corr)
  }
  
  # Convert per-gene Moran's I results into data.frames
  df_m676 <- do.call(rbind, moran_res676) %>% as.data.frame()
  df_m507 <- do.call(rbind, moran_res507) %>% as.data.frame()
  
  # Compute average SpaTM_S (corrected) Moran's I values per sample
  avg_m676 <- colMeans(df_m676, na.rm = TRUE)
  avg_m507 <- colMeans(df_m507, na.rm = TRUE)
  cat("Average Moran's I for 151676 (SpaTM_S):", avg_m676["SpaTM_S"], "\n")
  cat("Average Moran's I for 151507 (SpaTM_S):", avg_m507["SpaTM_S"], "\n")
  
  # Store these averages in the results matrix
  results[z, ] <- c(avg_m676["SpaTM_S"], avg_m507["SpaTM_S"])
  
  # Store detailed per-gene Moran's I values in master lists, keyed by gene set size
  all_moran676[[as.character(n_gene)]] <- df_m676
  all_moran507[[as.character(n_gene)]] <- df_m507
  
  models[[as.character(n_gene)]] <- list(spe_train, mlp)
  z <- z + 1
}
  
print(results)
write.csv(results, "./n_gene_fig2_benchmark_stm.csv")
saveRDS(all_moran676, "./all_moran676.rds")
saveRDS(all_moran507, "./all_moran507.rds")