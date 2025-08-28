set.seed(1234) 
library(SpaTM)
library(scran)
library(scater) 
library(tidyverse)   
library(Matrix) 
library(spatialLIBD)
library(torch)
library(DescTools)
library(mclust)
library(tictoc)
spe <- readRDS("~/Downloads/temp_brca.rds")


##Remove Mitochondrial genes
library(EnsDb.Hsapiens.v75)
gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")
spe <- spe[!rownames(spe) %in% gns$gene_name,]
## Remove Ribosomal Proteins
cur_genes <- genes(EnsDb.Hsapiens.v75)$symbol
cur_genes <- cur_genes[grepl(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",cur_genes)]
genes_remove <- genes(EnsDb.Hsapiens.v75, filter = ~ symbol %in% cur_genes)$gene_name
spe <- spe[!rownames(spe) %in% genes_remove,]





max_iter <- 100
lr <- 0.001

n_genes <- 5000
nbr_dist <- 3
spe_train <- spe
K <- 20
spe <- SpatialTopicExperiment(spe,K = K)

count_prop <- tabulate(counts(spe_train)@i + 1)/ncol(spe_train)
keep.genes <- rownames(spe_train)[count_prop >= 0.05 & count_prop <= 0.8]

##Get Top genes
spe_train <- computeGeometricFactors(spe_train)
spe_train <- logNormCounts(spe_train)
dec <- modelGeneVar(spe_train)
dec <- na.omit(dec)
dec <- dec[keep.genes,]
top_hvgs <- getTopHVGs(dec,n = n_genes)
spe_train <- spe_train[top_hvgs,]

#Remove cells & genes with no counts
spe_train <- spe_train[,colSums(counts(spe_train)) > 0]
spe_train <- spe_train[rowSums(counts(spe_train)) >0,]
top_hvgs <- rownames(spe_train)

spe_train <- SpatialTopicExperiment(
  spe_train,
  guided = FALSE,
  labels = NULL,
  K,
  hvg = NULL)


##Define Neighbors
spe_train$array_row <- spatialCoords(spe)[,1]
spe_train$array_col <- spatialCoords(spe)[,2]
nbr_list <- get_nbrs(spe_train,'sample_id','int_cell',NULL,dist = 3,loss_fun = 1)


#Train
D <- ncol(spe_train)


####Run RTM######
rtm_weights <- SpaTM:::train_RTM(counts(spe_train),
                         spe_train$int_cell,
                         rowData(spe_train)$gene_ints,
                         nbr_list,
                         alphaPrior(spe_train),
                         betaPrior(spe_train),
                         K,
                         ncol(spe_train),
                         ndk(spe_train),
                         nwk(spe_train),
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


#spe_out <- RTM(spe_train,K,nbr_list,1,6,100,lr,T,T,F,T)
spe_train <- buildTheta(spe_train)
spe_train <- buildPhi(spe_train)




anno_df <- read.csv("~/Downloads/BRCA_Data/metadata.tsv",header = T,sep = '\t')
rownames(anno_df) <- anno_df$ID
spe_test <- spe_train

spe_test$Anno <- anno_df[colnames(spe_test),]$fine_annot_type
##Build Adjacency Matrices 
metadata(spe_test)[['RTM_weights']] <- rtm_weights
reducedDim(spe_test,'ADJ') <- get_all_pred(spe_test,loss_fun = 1)
reducedDim(spe_test,'Theta') <- theta(spe_test)
num_clusts <- length(unique(spe_test$Anno))
spe_test <- runPCA(spe_test,dimred = 'ADJ')

### Generate Clusters
clust_input <- 'PCA'
library(bluster)
lei_clust <-  clust_optim(spe_test,num_clusts,200,X = clust_input, k = 10,alg = 'leiden')
lei_topic <- clust_optim(spe_test,num_clusts,200,X = 'Theta', k = 10,alg = 'leiden')
#repeat for louvain
lou_clust <-  clust_optim(spe_test,num_clusts,200,X = clust_input, k = 10,alg = 'louvain')
lou_topic <- clust_optim(spe_test,num_clusts,200,X = 'Theta', k = 10,alg = 'louvain')


spe_test$leiden <- lei_clust$clustering
spe_test$leiden <- lei_clust$clustering %>% as.character()
spe_test$lei_topic <- lei_topic$clustering %>% as.character()
spe_test$louvain <- lou_clust$clustering %>% as.character()
spe_test$lou_topic <- lou_topic$clustering %>% as.character()

spe_test$smooth <- rtm_smooth(spe_test,spe_test$leiden,nbr_list = nbr_list) 
spe_test$smooth <- rtm_smooth(spe_test,spe_test$smooth,k = 5)
spe_test$smoothlouvain <- rtm_smooth(spe_test,spe_test$louvain,nbr_list = nbr_list) 
spe_test$smoothlouvain <- rtm_smooth(spe_test,spe_test$smoothlouvain,k = 5)

spe_test$lei_topic <- rtm_smooth(spe_test,spe_test$lei_topic,nbr_list = nbr_list)
spe_test$lei_topic <- rtm_smooth(spe_test,spe_test$lei_topic,k = 5)
spe_test$lou_topic <- rtm_smooth(spe_test,spe_test$lou_topic,nbr_list = nbr_list)
spe_test$lou_topic <- rtm_smooth(spe_test,spe_test$lou_topic,k = 5)

spe_test$Topic <- apply(theta(spe_train),1,which.max)

spe_test$smoothTopic <- rtm_smooth(spe_test,spe_test$Topic,nbr_list = nbr_list) 
spe_test$smoothTopic <- rtm_smooth(spe_test,spe_test$smoothTopic,k = 5)

spe_test$pxl_col_in_fullres <- spe_test$imagecol
spe_test$pxl_row_in_fullres <- spe_test$imagerow
imgData(spe_test)$image_id <- 'lowres'
vis_clus(spe_test[,!is.na(spe_test$Anno)],clustervar = 'Anno')
vis_clus(spe_test,clustervar = 'leiden')
vis_clus(spe_test,clustervar = 'smooth')
vis_clus(spe_test,clustervar = 'smoothTopic')
vis_clus(spe_test,clustervar = 'lei_topic')
vis_clus(spe_test,clustervar = 'smoothlouvain')
vis_clus(spe_test,clustervar = 'lou_topic')
aricode::ARI(spe_test$Anno[!is.na(spe_test$Anno)],spe_test$smooth[!is.na(spe_test$Anno)])
aricode::ARI(spe_test$Anno[!is.na(spe_test$Anno)],spe_test$smoothTopic[!is.na(spe_test$Anno)])
aricode::ARI(spe_test$Anno[!is.na(spe_test$Anno)],spe_test$lei_topic[!is.na(spe_test$Anno)])
aricode::ARI(spe_test$Anno[!is.na(spe_test$Anno)],spe_test$lou_topic[!is.na(spe_test$Anno)])
aricode::ARI(spe_test$Anno[!is.na(spe_test$Anno)],spe_test$smoothlouvain[!is.na(spe_test$Anno)])


ggplot(as.data.frame(colData(spe_test)),aes(imagecol,-imagerow,color = Anno)) +
  geom_point()


viz_phi <- phi(spe_train)
#viz_phi <- t(apply(viz_phi,1,function(a){a/sum(a)}))
top10_k <- as.vector(sapply(1:ncol(viz_phi),function(K){
  top_10 <- order(viz_phi[,K],decreasing = T)[1:10]
  viz_phi[top_10,] <- 0
  top_10
}))

viz_phi <- phi(spe_train)[unique(top10_k),]



#col_scale <- colorRamp2(c(min(viz_phi), min(0.05,max(viz_phi))), c("white", "red"))
library(pheatmap)
library(circlize)
#col_scale <- colorRamp2(c(min(viz_phi), max(viz_phi)), c("white", "red"))
viz_phi[viz_phi > 0.05] <- 0.05
pheatmap(viz_phi,cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("white", "red"))(100),
         heatmap_legend_param = list(
           title = "Topic Weights",
           at = c(round(min(viz_phi,0)),
                  max(viz_phi))),
         cellwidth = 15, fontsize_row = 8)



vis_clus(spe_test,clustervar = 'Anno') + theme(legend.position = 'none')
for (i in 1:K){
  spe_test$Topic <- theta(spe_test)[,i]
  print(vis_gene(spe_test,geneid = 'Topic') + labs(title = (paste('Topic:',i))))
}


for (i in c(11,20,17,19)){
  spe_test$Topic <- theta(spe_test)[,i]
  print(vis_gene(spe_test,geneid = 'Topic') + labs(title = (paste('Topic:',i))))
}


#Topic 11
vis_gene(spe_test,geneid = 'TTLL12',minCount = -1)
vis_gene(spe_test,geneid = 'SREBF1',minCount = -1)
vis_gene(spe_test,geneid = 'THBS1',minCount = -1)
#Topic 20
vis_gene(spe_test,geneid = 'HMGB3',minCount = -1)
vis_gene(spe_test,geneid = 'S100A9',minCount = -1)
vis_gene(spe_test,geneid = 'CEBPD',minCount = -1)
#Topic 19
vis_gene(spe_test,geneid = 'IGHG2',minCount = -1)
vis_gene(spe_test,geneid = 'COMP',minCount = -1)
vis_gene(spe_test,geneid = 'CTGF',minCount = -1)
vis_gene(spe_test,geneid = 'MMP2',minCount = -1)
vis_gene(spe_test,geneid = 'SERPINF1',minCount = -1)
vis_gene(spe_test,geneid = 'ISLR',minCount = -1)
#Topic 17
vis_gene(spe_test,geneid = 'IGLC3',minCount = -1)
vis_gene(spe_test,geneid = 'JCHAIN',minCount = -1)
vis_gene(spe_test,geneid = 'CCL19',minCount = -1)
vis_gene(spe_test,geneid = 'TRBC2',minCount = -1)
vis_gene(spe_test,geneid = 'C1R',minCount = -1)
vis_gene(spe_test,geneid = 'COL6A3',minCount = -1)



library(fgsea)
library(GSEABase)
gmt_file <- "./Revisions/brca/h.all.v2025.1.Hs.symbols.gmt"
gene_set_collection <- getGmt(gmt_file)

# Convert to a named list
pathways <- lapply(gene_set_collection, geneIds)
names(pathways) <- sapply(gene_set_collection, function(x) x@setName)
results_list <- list()
gene_df <- as.data.frame(phi(spe_test))
gene_names <- rownames(gene_df)
for (factor_name in colnames(gene_df)) {
  cat("Running fgsea for:", factor_name, "\n")
  
  # Create named vector of gene scores (must be named by gene symbol)
  gene_scores <- gene_df[[factor_name]]
  names(gene_scores) <- gene_names
  
  # Sort scores in decreasing order
  gene_scores <- sort(gene_scores, decreasing = TRUE)
  
  # Run fgsea
  fgsea_res <- fgsea(pathways = pathways,
                     stats    = gene_scores,
                     minSize  = 15,
                     maxSize  = 500)
  
  # Save results
  fgsea_res$Topic <- factor_name
  results_list[[factor_name]] <- fgsea_res
}
res_df <- do.call('rbind',results_list)

