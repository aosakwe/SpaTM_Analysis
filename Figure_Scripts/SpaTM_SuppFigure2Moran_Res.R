#PLOT FOR FIGURE 3 

library(tidyverse)
library(reshape2)

#MORAN
files <- list.files("./Revisions/MORAN/")
files <- files[str_detect(files,"all_moran507")]
res <- data.frame(#Uncorrected = NA,
                  Moran = NA,
                  Method = NA,
                  NGENE = NA,
                  Gene = NA)
gene_vals <- c('500','1000','5000','10000','25000')
for (ngene in gene_vals){
  cur_files <- files[str_detect(files,paste0(ngene,'.rds'))]
  if (ngene == '5000'){
    cur_files <- cur_files[!str_detect(cur_files,paste0('25000','.rds'))]
  } else if (ngene == '1000'){
    cur_files <- cur_files[!str_detect(cur_files,paste0('10000','.rds'))]
  }
  cur <- readRDS(paste0("./Revisions/MORAN/",cur_files[str_detect(cur_files,'gtm')]))[[1]]
  cur2 <- readRDS(paste0("./Revisions/MORAN/",cur_files[!str_detect(cur_files,'gtm')]))[[1]]
  print(table(rownames(cur) == rownames(cur2)))
  cur$LDA <- cur2$SpaTM_S
  cur$NGENE <- ngene
  cur$Gene <- rownames(cur)
  cur <- melt(cur)
  cur <- cur[,c('value','variable','NGENE','Gene')]
  colnames(cur) <- colnames(res)
  rownames(cur) <- NULL
  res <- rbind(res,cur)
}
res <- res[-1,]
res$Method <- str_replace(res$Method,'SpaTM_S','SpaTM-S')

res$NGENE <- factor(res$NGENE,levels = c('500','1000','5000','10000','25000'))
ggplot(res,
       aes(NGENE,Moran,fill = Method)) +
  geom_boxplot() + 
  labs(y = "Moran's I",
       x = "# of Genes",
       title = 'Cross-Sample') +
  theme_minimal(base_size = 16) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5))



files <- list.files("./Revisions/MORAN/")
files <- files[str_detect(files,"all_moran676")]
res <- data.frame(#Uncorrected = NA,
  Moran = NA,
  Method = NA,
  NGENE = NA,
  Gene = NA)
gene_vals <- c('500','1000','5000','10000','25000')
for (ngene in gene_vals){
  cur_files <- files[str_detect(files,paste0(ngene,'.rds'))]
  if (ngene == '5000'){
    cur_files <- cur_files[!str_detect(cur_files,paste0('25000','.rds'))]
  } else if (ngene == '1000'){
    cur_files <- cur_files[!str_detect(cur_files,paste0('10000','.rds'))]
  }
  cur <- readRDS(paste0("./Revisions/MORAN/",cur_files[str_detect(cur_files,'gtm')]))[[1]]
  cur2 <- readRDS(paste0("./Revisions/MORAN/",cur_files[!str_detect(cur_files,'gtm')]))[[1]]
  print(table(rownames(cur) == rownames(cur2)))
  cur$LDA <- cur2$SpaTM_S
  cur$NGENE <- ngene
  cur$Gene <- rownames(cur)
  cur <- melt(cur)
  cur <- cur[,c('value','variable','NGENE','Gene')]
  colnames(cur) <- colnames(res)
  rownames(cur) <- NULL
  res <- rbind(res,cur)
}
res <- res[-1,]
res$Method <- str_replace(res$Method,'SpaTM_S','SpaTM-S')


res$NGENE <- factor(res$NGENE,levels = c('500','1000','5000','10000','25000'))
ggplot(res,
       aes(NGENE,Moran,fill = Method)) +
  geom_boxplot() + 
  labs(y = "Moran's I",
       x = "# of Genes",
       title = 'Within-Sample') +
  theme_minimal(base_size = 16) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5))



## Correlation
files <- list.files("./Revisions/MORAN/")
files <- files[str_detect(files,"151507")]
res <- list()
for (file in files){
  cur <- readRDS(paste0("./Revisions/MORAN/",file))
  res[[file]] <- cur
}
res <- do.call('cbind',res)
#colnames(res) <- str_replace(colnames(res),'n_gene_fig2_benchmark_stm_151507_','')
#colnames(res) <- str_replace(colnames(res),'.rds','')

res <- melt(res)
res$Method <- ifelse(str_detect(res$Var2,'gtm'),'LDA','SpaTM-S')
res$Var2 <- str_replace(res$Var2,'n_gene_fig2_benchmark_stm_151507_','')
res$Var2 <- str_replace(res$Var2,'n_gene_fig2_benchmark_gtm_151507_','')
res$Var2 <- str_replace(res$Var2,'.rds','')


res$Var2 <- factor(as.character(res$Var2),
                   levels = c('500','1000','5000','10000','25000'))
ggplot(res,
       aes(Var2,value,fill = Method)) + 
  geom_boxplot() +
  labs(x = "# of Genes",
       y = "Pearson's Correlation",
       title = 'Cross-Sample') +
  theme_minimal(base_size = 16) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5))


files <- list.files("./Revisions/MORAN/")
files <- files[str_detect(files,"151676")]
res <- list()
for (file in files){
  cur <- readRDS(paste0("./Revisions/MORAN/",file))
  res[[file]] <- cur
}
res <- do.call('cbind',res)
#colnames(res) <- str_replace(colnames(res),'n_gene_fig2_benchmark_stm_151507_','')
#colnames(res) <- str_replace(colnames(res),'.rds','')

res <- melt(res)
res$Method <- ifelse(str_detect(res$Var2,'gtm'),'LDA','SpaTM-S')
res$Var2 <- str_replace(res$Var2,'n_gene_fig2_benchmark_stm_151676_','')
res$Var2 <- str_replace(res$Var2,'n_gene_fig2_benchmark_gtm_151676_','')
res$Var2 <- str_replace(res$Var2,'.rds','')


res$Var2 <- factor(as.character(res$Var2),
                   levels = c('500','1000','5000','10000','25000'))
ggplot(res,
       aes(Var2,value,fill = Method)) + 
  geom_boxplot() +
  labs(x = "# of Genes",
       y = "Pearson's Correlation",
       title = 'Within-Sample') +
  theme_minimal(base_size = 16) +
  theme(axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5))
