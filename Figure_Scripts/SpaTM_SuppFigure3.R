##Benchmarking Script for SpaTM-R
set.seed(1234) 
library(SpaTM)
library(scran)
library(scater) 
library(tidyverse)   
library(Matrix) 
library(spatialLIBD)
library(mclust)
library(tictoc)
library(patchwork)
#source("./Figure_Scripts/alt_get_nbr.R")
#########################
get_nbrs_optimal <- function(spe, samples, cell_ids, group_by = NULL, 
                             dist_range = c(0.5, 5), target_median_nbrs = 10,
                             loss_fun = 1, max_iter = 20, tolerance = 0.5) {
  
  # Check array coordinates exist
  coldata <- colnames(colData(spe))
  if (!all(c('array_row','array_col') %in% coldata)) {
    stop('Error: array_row and array_col not found in spe metadata. This function assumes that spatial array coordinates are stored under these two names.')
  }
  
  # Compute distance matrix once
  dist_mat <- SpaTM:::get_dist_cpp(as.matrix(colData(spe)[,c('array_row','array_col')]))
  
  # Binary search for optimal distance threshold
  min_dist <- dist_range[1]
  max_dist <- dist_range[2]
  iter <- 0
  best_diff <- Inf
  best_nbr_list <- NULL
  best_dist <- NULL
  best_median <- NULL
  
  while (iter < max_iter && (max_dist - min_dist) > 0.01) {
    iter <- iter + 1
    curr_dist <- (min_dist + max_dist) / 2
    
    # Get neighbor list with current distance threshold
    nbr_list <- get_nbrs_with_dist(spe, samples, cell_ids, group_by, 
                                   dist_mat, curr_dist, loss_fun)
    
    # Calculate median number of neighbors
    nbr_counts <- sapply(nbr_list, function(x) {
      if (is.null(x) || nrow(x) == 0) return(0)
      sum(x[,2] == 1) # Count only positive neighbors
    })
    
    median_nbrs <- median(nbr_counts)
    diff_from_target <- abs(median_nbrs - target_median_nbrs)
    
    # Keep track of best result
    if (diff_from_target < best_diff) {
      best_diff <- diff_from_target
      best_nbr_list <- nbr_list
      best_dist <- curr_dist
      best_median <- median_nbrs
    }
    
    # If within tolerance, break early
    if (diff_from_target <= tolerance) {
      break
    }
    
    # Adjust search range
    if (median_nbrs < target_median_nbrs) {
      min_dist <- curr_dist
    } else {
      max_dist <- curr_dist
    }
    
    message(sprintf("Iter: %d, Dist: %.2f, Median nbrs: %.1f, Target: %d", 
                    iter, curr_dist, median_nbrs, target_median_nbrs))
  }
  
  if (is.null(best_nbr_list)) {
    warning("Failed to find optimal threshold. Using last computed threshold.")
    return(nbr_list)
  }
  
  message(sprintf("Optimization complete. Optimal distance: %.2f, Median neighbors: %.1f", 
                  best_dist, best_median))
  return(best_nbr_list)
}

get_nbrs_with_dist <- function(spe, samples, cell_ids, group_by, dist_mat, dist, loss_fun) {
  nbr_list <- list()
  
  if (!is.null(group_by)) {
    for (i in colData(spe)[,cell_ids]) {
      cur_sample <- colData(spe)[i,samples]
      cur_group <- colData(spe)[i,group_by]
      nbrs <- which(colData(spe)[,samples] == cur_sample & 
                      colData(spe)[,group_by] == cur_group &
                      dist_mat[i,] <= dist) - 1
      
      # Remove current example
      nbrs <- nbrs[nbrs != (i-1)]
      n_sam <- 50
      if (length(nbrs) > n_sam) {
        nbrs <- sample(nbrs, n_sam, 
                       prob = exp(-dist_mat[i,nbrs+1])/sum(exp(-dist_mat[i,nbrs+1])))
      }
      nbr_list[[i]] <- cbind(nbrs, rep(1, length(nbrs)))
      n <- length(nbrs)
      
      if (loss_fun == 1 & n > 0) {
        nbr_list[[i]] <- rbind(
          nbr_list[[i]],
          cbind(
            sample((1:ncol(spe))[-i], n, prob = 1/(1 + exp(-dist_mat[i,-i]))) - 1,
            rep(0, n)
          )
        )
      }
    }
  } else {
    for (i in colData(spe)[,cell_ids]) {
      cur_nbrs <- matrix(0, nrow = 0, ncol = 2)
      cur_sample <- colData(spe)[i,samples]
      nbrs <- which(dist_mat[i,] <= dist & colData(spe)[,samples] == cur_sample) - 1
      
      # Remove current example
      nbrs <- nbrs[nbrs != (i-1)]
      n <- length(nbrs)
      nbr_list[[i]] <- rbind(cur_nbrs, cbind(nbrs, rep(1, n)))
      
      if (loss_fun != 0 & n > 0) {
        if (loss_fun == 1) {
          nbr_list[[i]] <- rbind(
            nbr_list[[i]],
            cbind(
              sample((1:ncol(spe))[-i], n, prob = 1/(1 + exp(-dist_mat[i,-i]))) - 1,
              rep(0, n)
            )
          )
        } else if (loss_fun == 2) {
          # Euclidean Distance Prediction
          nbr_list[[i]] <- rbind(
            nbr_list[[i]],
            cbind(sample((1:ncol(spe))[-i], n) - 1, rep(0, n))
          )
          nbr_list[[i]] <- cbind(
            nbr_list[[i]],
            log(dist_mat[i, nbr_list[[i]][,1] + 1])
          )
        }
      } else if (n == 0) {
        nbr_list[[i]] <- rbind(cbind(cur_nbrs, rep(0, n)), cbind(nbrs, rep(1, n), rep(0, n)))
      }
    }
  }
  
  if (do.call('sum', lapply(nbr_list, nrow)) == 0) {
    message('Warning. Your current parametrization returned 0 neighbors for all samples. Consider a different distance threshold or grouping by a covariate.')
  }
  
  return(nbr_list)
}
print('here')
get_nbrs_optimal <- function(spe, samples, cell_ids, group_by = NULL, 
                             k = 10, loss_fun = 1) {
  # Check array coordinates exist
  coldata <- colnames(colData(spe))
  if (!all(c('array_row','array_col') %in% coldata)) {
    stop('Error: array_row and array_col not found in spe metadata. This function assumes that spatial array coordinates are stored under these two names.')
  }
  
  # Compute distance matrix once
  dist_mat <- SpaTM:::get_dist_cpp(as.matrix(colData(spe)[, c('array_row','array_col')]))
  
  nbr_list <- list()
  n_cells <- ncol(spe)
  
  if (!is.null(group_by)) {
    for (i in seq_len(n_cells)) {
      cur_sample <- colData(spe)[i, samples]
      cur_group <- colData(spe)[i, group_by]
      # Select indices within the same sample and group
      valid_inds <- which(colData(spe)[, samples] == cur_sample & colData(spe)[, group_by] == cur_group)
      valid_inds <- setdiff(valid_inds, i)  # remove self
      if (length(valid_inds) == 0) {
        nbr_list[[i]] <- matrix(numeric(0), ncol = 2)
        next
      }
      # Order based on distance
      dists <- dist_mat[i, valid_inds]
      nn_inds <- valid_inds[order(dists)][1:min(k, length(valid_inds))]
      pos_nbrs <- cbind(nn_inds - 1, rep(1, length(nn_inds)))
      
      # Optionally add negative neighbors if loss_fun == 1
      if (loss_fun == 1) {
        all_inds <- setdiff(seq_len(n_cells), i)
        neg_pool <- setdiff(all_inds, nn_inds)
        neg_sample_size <- min(length(nn_inds), length(neg_pool))
        if (neg_sample_size > 0) {
          neg_inds <- sample(neg_pool, neg_sample_size)
          neg_nbrs <- cbind(neg_inds - 1, rep(0, neg_sample_size))
          nbr_list[[i]] <- rbind(pos_nbrs, neg_nbrs)
        } else {
          nbr_list[[i]] <- pos_nbrs
        }
      } else {
        nbr_list[[i]] <- pos_nbrs
      }
    }
  } else {
    for (i in seq_len(n_cells)) {
      valid_inds <- setdiff(seq_len(n_cells), i)
      dists <- dist_mat[i, valid_inds]
      nn_inds <- valid_inds[order(dists)][1:min(k, length(valid_inds))]
      pos_nbrs <- cbind(nn_inds - 1, rep(1, length(nn_inds)))
      
      # Optionally add negative neighbors if loss_fun == 1
      if (loss_fun == 1) {
        neg_pool <- setdiff(valid_inds, nn_inds)
        neg_sample_size <- min(length(nn_inds), length(neg_pool))
        if (neg_sample_size > 0) {
          neg_inds <- sample(neg_pool, neg_sample_size)
          neg_nbrs <- cbind(neg_inds - 1, rep(0, neg_sample_size))
          nbr_list[[i]] <- rbind(pos_nbrs, neg_nbrs)
        } else {
          nbr_list[[i]] <- pos_nbrs
        }
      } else {
        nbr_list[[i]] <- pos_nbrs
      }
    }
  }
  
  if (do.call('sum', lapply(nbr_list, nrow)) == 0) {
    message('Warning. The parametrization returned 0 neighbors for all samples.')
  }
  
  return(nbr_list)
}
######################
all_sce <- list.files("../data/SpatialClusterBenchmark/")
all_sce <- all_sce[str_detect(all_sce,'rds')]
max_iter <- 100
lr <- 1e-5
for(file in all_sce){
  spe <- readRDS(paste0("../data/SpatialClusterBenchmark/",file))
  spe$sample_id <- str_replace(file,'.rds','')
  K <- length(unique(spe$ground_truth))
  if( nrow(spe) < 3000){
    spe <- SingleCellTopicExperiment(spe,
                                          K = K,
                                          hvg = NULL)
  } else {
    spe <- SingleCellTopicExperiment(spe,
                                          K = K,
                                          hvg = 3000)
  }
  ##########################
  #nbr_list <- get_nbrs(spe,'sample_id','int_cell',NULL,dist = 50,loss_fun = 1)
  spe <- spe[,colSums(counts(spe)) > 0]
  spe <- SingleCellTopicExperiment(spe,
                                   K = K,
                                   hvg = NULL)
  nbr_list <- get_nbrs_optimal(spe,
                               'sample_id',
                               'int_cell',
                               k = 10,
                               loss_fun = 1)
  #Train
  D <- ncol(spe)
  
  
  ####Run RTM######
  #tic()
  rtm_weights <- SpaTM:::train_RTM(counts(spe),
                                   spe$int_cell,
                                   rowData(spe)$gene_ints,
                                   nbr_list,
                                   alphaPrior(spe),
                                   betaPrior(spe),
                                   K,
                                   ncol(spe),
                                   ndk(spe),
                                   nwk(spe),
                                   6,
                                   max_iter,
                                   verbal = TRUE,
                                   zero_gamma = FALSE,
                                   rand_gamma = TRUE,#FALSE,
                                   thresh = 0.001,
                                   lr =lr,
                                   rho = 50000,
                                   loss_fun = 1,
                                   m_update = TRUE)
  #toc()
  
  #spe_out <- RTM(spe,K,nbr_list,1,6,100,lr,T,T,F,T)
  spe <- buildTheta(spe)
  spe <- buildPhi(spe)
  
  
  
  
  ##Build Adjacency Matrices 
  metadata(spe)[['RTM_weights']] <- rtm_weights
  reducedDim(spe,'ADJ') <- get_all_pred(spe,loss_fun = 1)
  reducedDim(spe,'Theta') <- theta(spe)
  num_clusts <- length(unique(spe$ground_truth))
  spe <- runPCA(spe,dimred = 'ADJ')
  
  ### Generate Clusters
  clust_input <- 'PCA'
  library(bluster)
  lei_clust <-  clust_optim(spe,num_clusts,200,X = clust_input, k = 10,alg = 'louvain')
  lei_topic <- clust_optim(spe,num_clusts,200,X = 'Theta', k = 10,alg = 'louvain')
  
  spe$leiden <- lei_clust$clustering
  spe$leiden <- lei_clust$clustering %>% as.character()
  spe$lei_topic <- lei_topic$clustering %>% as.character()
  
  spe$smooth <- rtm_smooth(spe,spe$leiden,nbr_list = nbr_list) 
  spe$smooth <- rtm_smooth(spe,spe$smooth,k = 5)
  spe$lei_topic <- rtm_smooth(spe,spe$lei_topic,nbr_list = nbr_list)
  spe$lei_topic <- rtm_smooth(spe,spe$lei_topic,k = 5)
  
  spe$Topic <- apply(theta(spe),1,which.max)
  
  spe$smoothTopic <- rtm_smooth(spe,spe$Topic,nbr_list = nbr_list) 
  spe$smoothTopic <- rtm_smooth(spe,spe$smoothTopic,k = 5)
  saveRDS(spe,paste0("../data/SpatialClusterBenchmark/spatm_",file))
  #############################
  
  
  p1 <- colData(spe) %>%
    as.data.frame() %>%
    ggplot(aes(array_row,array_col,color = smoothTopic)) +
    geom_point() + labs(title = file)
  p2 <- colData(spe) %>%
    as.data.frame() %>%
    ggplot(aes(array_row,array_col,color = smooth)) +
    geom_point() + labs(title = file)
  p3 <- colData(spe) %>%
    as.data.frame() %>%
    ggplot(aes(array_row,array_col,color = lei_topic)) +
    geom_point() + labs(title = file)
  p4 <- colData(spe) %>%
    as.data.frame() %>%
    ggplot(aes(array_row,array_col,color = ground_truth)) +
    geom_point() + labs(title = file)
  print(wrap_plots(list(p1,p2,p3,p4),ncol = 2))
  ggsave(paste(file,"RTM.png",sep = "_"),width = 20,height = 20)
}

