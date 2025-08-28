##SlIDESEQ eVALUATION
set.seed(1234) 
library(SpaTM)
library(scran)
library(scater) 
library(tidyverse)   
library(Matrix) 
library(spatialLIBD) 
library(torch)


sce <- readRDS("./slide_sce.rds")
#spe <- spe_train
spe <- readRDS("./Revisions/SlideTags/SpaTM_SlideSeq_Train_03JulyBalAlt_1000_8.rds")
spe <- spe$spe_train
rownames(phi(spe)) <- rowData(spe)$gene_name
mlp <- torch_load("./Revisions/SlideTags/SpaTM_SlideSeq_Train_03JulyBalAlt1000_8.torch")


keep.genes <- intersect(rownames(sce),rownames(phi(spe)))
sce <- sce[keep.genes,]
phi(spe) <- phi(spe)[keep.genes,]


sce <- SingleCellTopicExperiment(sce,
                                 K  = ncol(phi(spe)))
sce <- inferTopics(sce,
                   num_thread = 6,
                   maxiter = 150,
                   verbal = T,
                   phi = phi(spe))
sce <- buildTheta(sce)
sce$Pred <- levels(spe$spatialLIBD)[pred_mlp(mlp,theta(sce))]
sce$layer <- sce$Pred
table(sce$cell_type,sce$Pred)
df <- read.csv("./data/SlideSeqDLPFC/humancortex_spatial.csv")
df <- df[-1,]
colnames(sce) <- str_replace(colnames(sce),'\\.','-')
table(colnames(sce) == df$NAME)

sce$X <- as.numeric(df$X)
sce$Y <- as.numeric(df$Y)

reducedDim(sce,"SPATIAL") <- colData(sce)[,c('X','Y')]
sce$Layer <- factor(sce$layer,
                    levels = c('L1','L2','L3','L4','L5','L6','WM'))



layer_colors <- c(
  'L1' = "#b2df8a",
  'L2' = "#e41a1c",
  'L3' = "#377eb8",
  'L4' = "#4daf4a",
  'L5' = "#ff7f00",
  'L6' = "gold",
  'WM' = "#a65628"
)

plotReducedDim(sce,'SPATIAL',
               colour_by = 'Layer') +
  scale_color_manual(values = layer_colors) +
  geom_point(alpha = 0.001)




df <- cbind(reducedDim(sce,'SPATIAL'),colData(sce))
ggplot(as.data.frame(df),aes(X,Y,color = Layer)) +
  geom_point() +
  facet_wrap( ~ Layer)


library(FNN) 

# Extract spatial coordinates from reducedDims
spatial_coords <- reducedDim(sce, "SPATIAL")

# Extract Layer labels
layers <- colData(sce)$Layer

# Number of neighbors (choose based on your resolution, e.g., 10)
k <- 10

# Get indices of K nearest neighbors for each cell
knn_result <- get.knn(spatial_coords, k = k)
knn_indices <- knn_result$nn.index

# Function to get the majority layer label among neighbors
get_majority <- function(indices) {
  neighbor_labels <- layers[indices]
  # Return most frequent label (mode)
  names(sort(table(neighbor_labels), decreasing = TRUE))[1]
}

# Apply smoothing across all cells
smoothed_layer <- apply(knn_indices, 1, get_majority)

# Add new column to colData
colData(sce)$Knn <- smoothed_layer


plotReducedDim(sce,'SPATIAL',
               colour_by = 'Layer') +
  scale_color_manual(values = layer_colors) +
  geom_point(alpha = 0.001)

plotReducedDim(sce,'SPATIAL',
               colour_by = 'Layer') +
  geom_point(alpha = 0.001)

plotReducedDim(sce,'SPATIAL',
               colour_by = 'Knn') +
  scale_color_manual(values = layer_colors) +
  geom_point(alpha = 0.001)

plotReducedDim(sce,'SPATIAL',
               colour_by = 'cell_type') +
  geom_point(alpha = 0.001)


