rm(list = ls())
###文件输入######
dat <- data.table::fread("Disease_Matrix.csv",data.table = F)
rownames(dat) <- dat$V1
dat <- dat[,-1]
# genes <- readxl::read_xlsx("hubGenes.xlsx")
genes <- data.table::fread("ModuleGenes.csv",data.table = F)

dat_filter <- as.data.frame(dat[genes$x,])
dat_filter = as.matrix(dat_filter)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")

###Consensus Cluster Analysis######
library(ConsensusClusterPlus)
cluster <- ConsensusClusterPlus(dat_filter,
                             maxK=9,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title='Consensus_Cluster',
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")
# clusterAlg="km",distance="euclidean"/"pam"+"pearson"/"pam"+"spearman"
#确认分组情况
library(dplyr)
label <- cluster[[2]][["consensusClass"]] %>% as.data.frame()
colnames(label) <- "cluster"
label$sample <- rownames(label)
label <- label[order(label$cluster),]
dat_Cluster <- dat[genes$x,rownames(label)]
dat_hub_boxplot <- t(rbind(label$cluster,dat_Cluster[genes$x,]))
colnames(dat_hub_boxplot)[1] <- "Cluster"
dat_hub_boxplot <- as.data.frame(dat_hub_boxplot)
dat_hub_boxplot$Cluster <- gsub("1", "Cluster1", dat_hub_boxplot$Cluster, fixed = TRUE)#把分组中的具体分组名称重新命名
dat_hub_boxplot$Cluster <- gsub("2", "Cluster2", dat_hub_boxplot$Cluster, fixed = TRUE)
# dat_hub_boxplot$Cluster <- gsub("3", "Cluster3", dat_hub_boxplot$Cluster, fixed = TRUE)
write.csv(dat_hub_boxplot,"6-Boxplot.csv", row.names = F)

dat_heatmap <- dat_Cluster[,rownames(label)]
dat_heatmap <- rbind(dat_hub_boxplot$Cluster,colnames(dat_heatmap),dat_heatmap)
rownames(dat_heatmap)[1] <- "#Group"
rownames(dat_heatmap)[2] <- "id"
write.table(dat_heatmap,"5-Heatmap.csv",col.names = F,quote = F ,sep = ",")
###Consensus Cluster Analysis可视化######
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)

length(which(label$cluster %in% "1")) #39
length(which(label$cluster %in% "2")) #55
# length(which(label$cluster %in% "3")) #16
# length(which(label$cluster %in% "4")) #16
c1 <- label$sample[which(label$cluster %in% "1")]
c2 <- label$sample[which(label$cluster %in% "2")]
# c3 <- label$sample[which(label$cluster %in% "3")]
# c4 <- label$sample[which(label$cluster %in% "4")]
# write.table(c1,"1-Cluster1.txt",row.names = F,col.names = F,quote = F)
# write.table(c2,"2-Cluster2.txt",row.names = F,col.names = F,quote = F)
# write.table(c3,"3-Cluster3.txt",row.names = F,col.names = F,quote = F)
# write.table(c4,"4-Cluster4.txt",row.names = F,col.names = F,quote = F)
annotation_col = data.frame(Type = factor(c(rep('Cluster1',39),rep('Cluster2',55))))
rownames(annotation_col) = colnames(dat_Cluster)

###2D PCA图######
# library(FactoMineR)
# library(factoextra)
# iris.pca <- PCA(t(dat_Cluster),graph = F)
# cla <- c(rep('Cluster1',39),rep('Cluster2',55))
# pdf("4-Cluster_PCA.pdf",width =8.1,height =5)
# ind.p <- fviz_pca_ind(iris.pca,geom.ind = "point", col.ind = cla,
#                     palette = c("#f0a274", "#2486b9"),
#                     legend.title = "Cluster",addEllipses = TRUE)
# ggpubr::ggpar(ind.p,title = "Consensus Cluster PCA",ggtheme = theme_bw())
# dev.off()

###3D PCA图######
library(Rtsne)
cla <- c(rep('Cluster1',39),rep('Cluster2',55))
tSNE_res <- Rtsne(t(dat_Cluster),
                  dims=3,
                  perplexity=5,
                  verbose=F,
                  max_iter=500,
                  check_duplicates=F)
tsne <- data.frame(tSNE1 = tSNE_res[["Y"]][,1],
                   tSNE2 = tSNE_res[["Y"]][,2],
                   tSNE3 = tSNE_res[["Y"]][,3],
                   cluster = cla)

tsne$colors <- ifelse(tsne$cluster %in% "Cluster1","#f0a274","#2486b9")
colnames(tsne) <- c("x","y","z","cluster","colors")
tsne$cluster <- as.factor(tsne$cluster)
pdf("4-Cluster_PCA.pdf",width = 8.1,height = 6)
library(scatterplot3d)
p3 <- scatterplot3d(tsne[,1:3],color = tsne$colors,main="Consensus Cluster PCA",pch = 21,bg = tsne$colors)
legend("bottom",col = "black", legend = levels(tsne$cluster),pt.bg = c("#f0a274","#2486b9"), pch = 21,
       inset = -0.2, xpd = TRUE, horiz = TRUE)
dev.off()

write.csv(dat_hub_boxplot[,1,drop = F],"Cluster_Group.csv",row.names = T)
write.csv(dat_hub_boxplot[,1,drop = F],"../14-Cluster_Immune/Cluster_Group.csv",row.names = T)
###END#####