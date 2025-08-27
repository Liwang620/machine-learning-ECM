rm(list=ls())
library(GEOquery)

###GSE75010######
##Download Datasets
aset <- getGEO("GSE75010",destdir = ".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)  ## 平台文件 
# exprs()提取matrix；pData提取group info
GSE75010 <- exprs(aset$GSE75010_series_matrix.txt.gz)
GSE75010 <- as.data.frame(GSE75010)
PD75010 <- pData(aset$GSE75010_series_matrix.txt.gz)

##GSE Platform
library(dplyr)
aset$GSE75010_series_matrix.txt.gz$platform_id #GPL6244
GPL6244 <- getGEO("GPL6244")
anno = GPL6244@dataTable@table

##Group Information
#把数据集GSE75010矩阵和分组情况保存下来
PD75010 <- PD75010[,c(2,10)]
colnames(PD75010)[2] <- 'Group'
table(PD75010$Group)

PD75010$Group <- gsub("diagnosis: non-PE", 
                      "Control", PD75010$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD75010$Group <- gsub("diagnosis: PE", 
                      "PE", PD75010$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD75010 <- PD75010[order(PD75010$Group,decreasing = F),]

#把pd里和dat里取交集 然后糅合到一起 删除dat多余的没pd那些样本的数据（因为pd只有一些数据，dat有很多。）
GSE75010 <- GSE75010[,match(PD75010$geo_accession,colnames(GSE75010))] 
colnames(anno)

library(magrittr)
anno$GeneSymbol <- strsplit(anno$gene_assignment, split = " /// ") %>% lapply(function(x){
  # x = strsplit(dat$gene_assignment[70750:70753], split = " /// ")[[1]]
  x1 <- x
  if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
  if(length(x1) == 0) x1 <- x[1]
  vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
    unlist() %>% unique()
  vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
  vec_gene <- paste0(vec_gene, collapse = " /// ")
  return(vec_gene)
}) %>% unlist()
anno$GeneSymbol %>% table() %>% sort(decreasing = T) %>% .[1:10]

GPLanno <- anno[,c("ID", "GeneSymbol")]
GSE75010$ID = rownames(GSE75010) #新建了一个ID用于当做新的一个列 

GSE75010 <- merge(GPLanno,GSE75010, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
GSE75010 <- na.omit(GSE75010)
GSE75010 <- aggregate(x=GSE75010,by = GSE75010$`GeneSymbol` %>% list(),FUN = mean)#多个探针对应一个基因取平均
rownames(GSE75010) = GSE75010$Group.1 #把Symbol变成行名
GSE75010 <- GSE75010[-1,-(1:3)]
GSE75010 <- GSE75010[!grepl("///", row.names(GSE75010)),]
range(GSE75010)
# GSE75010 <- log2(GSE75010+1)

###GSE30186######
##Download Datasets
bset <- getGEO("GSE30186",destdir = ".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)  ## 平台文件 
# exprs()提取matrix；pData提取group info
GSE30186 <- exprs(bset$GSE30186_series_matrix.txt.gz)
GSE30186 <- as.data.frame(GSE30186)
PD30186 <- pData(bset$GSE30186_series_matrix.txt.gz)

##GSE Platform
bset$GSE30186_series_matrix.txt.gz$platform_id #GPL10558
GPL10558 <- getGEO("GPL10558")
bnno = GPL10558@dataTable@table

##Group Information
#把数据集GSE30186矩阵和分组情况保存下来
PD30186 <- PD30186[,c(2,8)]
colnames(PD30186)[2] <- 'Group'
table(PD30186$Group)

PD30186$Group <- gsub("normal", 
                      "Control", PD30186$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD30186$Group <- gsub("preeclampsia", 
                      "PE", PD30186$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD30186 <- PD30186[order(PD30186$Group,decreasing = F),]
# PD30186 <- PD30186[-(9:15),]

#把pd里和dat里取交集 然后糅合到一起 删除dat多余的没pd那些样本的数据（因为pd只有一些数据，dat有很多。）
GSE30186 <- GSE30186[,match(PD30186$geo_accession,colnames(GSE30186))] 
colnames(bnno)

GPLbnno <- bnno[,c("ID", "Symbol")] 
GSE30186$ID = rownames(GSE30186) #新建了一个ID用于当做新的一个列 

GSE30186 <- merge(GPLbnno,GSE30186, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
GSE30186 <- na.omit(GSE30186)
GSE30186 <- aggregate(x=GSE30186,by = GSE30186$Symbol %>% list(),FUN = mean)#多个探针对应一个基因取平均
rownames(GSE30186) = GSE30186$Group.1 #把Symbol变成行名
GSE30186 <- GSE30186[-1,-(1:3)]
GSE30186 <- GSE30186[!grepl("///", row.names(GSE30186)),]
range(GSE30186)
GSE30186 = GSE30186-apply(GSE30186,1,min) #包含负值
# GSE30186 <- GSE30186-min(GSE30186)
GSE30186 <- log2(GSE30186+1)

###GSE24129######
##Download Datasets
cset <- getGEO("GSE24129",destdir = ".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)  ## 平台文件 
# exprs()提取matrix；pData提取group info
GSE24129 <- exprs(cset$GSE24129_series_matrix.txt.gz)
GSE24129 <- as.data.frame(GSE24129)
PD24129 <- pData(cset$GSE24129_series_matrix.txt.gz)

##GSE Platform
library(dplyr)
cset$GSE24129_series_matrix.txt.gz$platform_id #GPL6244
# GPL6244 <- getGEO("GPL6244")
cnno = GPL6244@dataTable@table

##Group Information
#把数据集GSE24129矩阵和分组情况保存下来
PD24129 <- PD24129[,c(2,11)]
colnames(PD24129)[2] <- 'Group'
table(PD24129$Group)

PD24129$Group <- gsub("disease state: Normotensive control pregnancy", 
                      "Control", PD24129$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD24129$Group <- gsub("disease state: Pre-eclampsia", 
                      "PE", PD24129$Group, fixed = TRUE)#把分组中的具体分组名称重新命名
PD24129 <- PD24129[order(PD24129$Group,decreasing = F),]
# grep("Fetal growth restriction",PD24129$Group)
PD24129 <- PD24129[-c(grep("Fetal growth restriction",PD24129$Group)),]

#把pd里和dat里取交集 然后糅合到一起 删除dat多余的没pd那些样本的数据（因为pd只有一些数据，dat有很多。）
GSE24129 <- GSE24129[,match(PD24129$geo_accession,colnames(GSE24129))] 
colnames(cnno)

# library(magrittr)
# anno$GeneSymbol <- strsplit(anno$gene_assignment, split = " /// ") %>% lapply(function(x){
#   # x = strsplit(dat$gene_assignment[70750:70753], split = " /// ")[[1]]
#   x1 <- x
#   if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
#   if(length(x1) == 0) x1 <- x[1]
#   vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
#     unlist() %>% unique()
#   vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
#   vec_gene <- paste0(vec_gene, collapse = " /// ")
#   return(vec_gene)
# }) %>% unlist()
# anno$GeneSymbol %>% table() %>% sort(decreasing = T) %>% .[1:10]

GPLcnno <- anno[,c("ID", "GeneSymbol")]
GSE24129$ID = rownames(GSE24129) #新建了一个ID用于当做新的一个列 

GSE24129 <- merge(GPLcnno,GSE24129, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
GSE24129 <- na.omit(GSE24129)
GSE24129 <- aggregate(x=GSE24129,by = GSE24129$`GeneSymbol` %>% list(),FUN = mean)#多个探针对应一个基因取平均
rownames(GSE24129) = GSE24129$Group.1 #把Symbol变成行名
GSE24129 <- GSE24129[-1,-(1:3)]
GSE24129 <- GSE24129[!grepl("///", row.names(GSE24129)),]
range(GSE24129)
# GSE24129 <- log2(GSE24129+1)

###Combined Datasets######
##寻找三个个数据集中都有的基因
genes <- Reduce(intersect,list(rownames(GSE75010),rownames(GSE30186),rownames(GSE24129)))
# genes <- intersect(rownames(GSE75010),rownames(GSE30186))
GSE75010 <- GSE75010[genes,]
GSE30186 <- GSE30186[genes,]
GSE24129 <- GSE24129[genes,]

#！！！注意，此处combined的数据集顺序要一致
combined_gse <- cbind(GSE75010,GSE30186,GSE24129)
combined_pd <- rbind(PD75010,PD30186,PD24129)

###画校准前后的Boxplot和PCA######
##Set Theme
library(magrittr)
library(ggplot2)
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # 图像
        plot.title = element_text(size = 7)) + 
  theme(panel.background = element_blank(), # 面板
        panel.grid = element_blank(),
        panel.borde = element_rect(fill = NA, size = 0.75 * 0.47)) + # 添加外框
  theme(axis.line = element_line(size = 0.75 * 0.47), # 坐标轴
        axis.text = element_text(size = 6, color = "black"), # 坐标文字为黑色
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.75 * 0.47)) +
  theme(legend.key = element_rect(fill = "white"),
        legend.key.size = unit(c(0.3, 0.3), "cm"), # 图注
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())

##画矫正前的boxplot图
#转长数据
dat_before <- cbind(rownames(combined_gse),combined_gse)
rownames(dat_before) <- NULL
colnames(dat_before)[1] <- "Gene"

sample_group <- combined_pd
rownames(sample_group) <- NULL
table(sample_group$sample)#查看样本数量
#！！！注意，此处赋值的顺序要和combined_pd的顺序要一直
sample_group$Group <- c(rep("GSE75010",157),rep("GSE30186",12),rep("GSE24129",16))

dat_before_long <- dat_before %>% 
  tibble::column_to_rownames("Gene") %>% t() %>% 
  data.frame() %>% 
  dplyr::mutate(sample = rownames(.)) %>% 
  merge(sample_group, by.x = "sample", by.y = "geo_accession") %>% 
  tidyr::gather(key = geneid, value, -c(sample, Group))

#Before Boxplot绘图
boxplot_before <- 
  ggplot(dat_before_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = Group), lwd = 0.1 *0.47, outlier.shape = NA) + # 以分组填充颜色，不要离群点
  scale_fill_manual(values = c("GSE75010" = "#F0A274", "GSE30186" = "#83a78d", "GSE24129" = "#eea2a4")) +
  labs(title = "Before Normalization") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank(), # 隐藏y轴标签
        axis.text.x = element_blank(), # 隐藏x轴内容
        axis.ticks.x = element_blank()) + # 隐藏x轴刻度线
  scale_y_continuous(limits = c(0,15),
                    # breaks = c(seq(150, 450,by = 50)), 
                     expand = c(0,0)) +
  coord_cartesian(clip = "off")# 解决上右框线看起来浅的问题
boxplot_before
ggsave("1-Boxpot_Before.pdf", plot = boxplot_before, units = "cm",width = 8.1,height = 5) 

#BiocManager::install("sva")
#数据矫正
library(sva)
#library(tidyverse)
library(limma)
#  "1"代表GSE75010,"2"代表GSE30186,"3"代表GSE24129
#！！！注意，此处赋值顺序与combined_gse顺序一致
batch <- c(rep("1",157),rep("2",12),rep("3",16))
#去除批次效应
adjusted_gse <- ComBat(combined_gse, batch=batch) 
#标准化
adjusted_gse <- normalizeBetweenArrays(adjusted_gse)

#画矫正后的boxplot图
#转长数据
adjusted_gse <- as.data.frame(adjusted_gse)
dat_after <- cbind(rownames(adjusted_gse),adjusted_gse)
rownames(dat_after) <- NULL
colnames(dat_after)[1] <- "Gene"
  
dat_after_long <- dat_after %>% 
  tibble::column_to_rownames("Gene") %>% t() %>% 
  data.frame() %>% 
  dplyr::mutate(sample = rownames(.)) %>% 
  merge(sample_group, by.x = "sample", by.y = "geo_accession") %>% 
  tidyr::gather(key = geneid, value, -c(sample, Group))

#After Boxplot绘图
boxplot_after <- 
  ggplot(dat_after_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = Group), lwd = 0.1 *0.47, outlier.shape = NA) + # 以分组填充颜色，不要离群点
  scale_fill_manual(values = c("GSE75010" = "#F0A274", "GSE30186" = "#83a78d", "GSE24129" = "#eea2a4")) +
  labs(title = "After Normalization") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank(), # 隐藏y轴标签
        axis.text.x = element_blank(), # 隐藏x轴内容
        axis.ticks.x = element_blank()) + # 隐藏x轴刻度线
  scale_y_continuous(limits = c(0,15),
                     # breaks = c(seq(150, 450,by = 50)), 
                     expand = c(0,0)) +
  coord_cartesian(clip = "off")
boxplot_after
ggsave("2-Boxplot_After.pdf", plot = boxplot_after, units = "cm",width = 8.1,height = 5) 


library(FactoMineR)
library(factoextra)
#绘制矫正前的PCA图
matrix_before <- as.matrix(t(combined_gse))
cla <- c(rep("GSE75010",157),rep("GSE30186",12),rep("GSE24129",16))
iris.pca_before <- PCA(matrix_before,graph = F,scale.unit = TRUE)

pdf("3-PCA_Before.pdf",width = 6,height = 4.5)
ind.p_before <- fviz_pca_ind(iris.pca_before,geom.ind = "point", col.ind = cla,
                    palette = c("GSE75010" = "#F0A274", "GSE30186" = "#83a78d", "GSE24129" = "#eea2a4"),
                    addEllipses = TRUE,
                    legend.title = "Group" )
ggpubr::ggpar(ind.p_before,title = "Before Normalization",ggtheme = theme_bw())
dev.off()
# #矫正前的3D PCA图
# library(Rtsne)
# tSNE_res <- Rtsne(matrix_before,
#                   dims=3,
#                   perplexity=5,
#                   verbose=F,
#                   max_iter=500,
#                   check_duplicates=F)
# tsne <- data.frame(tSNE1 = tSNE_res[["Y"]][,1],
#                    tSNE2 = tSNE_res[["Y"]][,2],
#                    tSNE3 = tSNE_res[["Y"]][,3],
#                    cluster = cla)
# 
# tsne$colors <- ifelse(tsne$cluster %in% "GSE75010","#1D7E8F",ifelse(tsne$cluster %in% "GSE30186","#F0A274","#B8DADC"))
# colnames(tsne) <- c("x","y","z","cluster","colors")
# tsne$cluster <- as.factor(tsne$cluster)
# # pdf("3-PCA_Before.pdf",width = 8.1,height = 6)
# library(scatterplot3d)
# pca_before <- scatterplot3d(tsne[,1:3],color = tsne$colors,main="Before Normalization",pch = 21,bg = tsne$colors)
# legend("bottom",col = "black", legend = levels(tsne$cluster),
#        pt.bg = c("#B8DADC","#1D7E8F","#F0A274"), 
#        pch = 21, inset = -0.2, xpd = TRUE, horiz = TRUE)
# dev.off()

#绘制矫正后的PCA图
matrix_after <- as.matrix(t(adjusted_gse))
iris.pca_after <- PCA(matrix_after,graph =F)
pdf("4-PCA_After.pdf",width = 6,height = 4.5)
ind.p_after <- fviz_pca_ind(iris.pca_after,geom.ind = "point", col.ind = cla,
                    palette = c("GSE75010" = "#F0A274", "GSE30186" = "#83a78d", "GSE24129" = "#eea2a4"),
                    addEllipses = TRUE,
                    legend.title = "Group")
ggpubr::ggpar(ind.p_after,title = "After Normalization",ggtheme = theme_bw())
dev.off()

# #矫正后的3D PCA图
# tSNE_res <- Rtsne(matrix_after,
#                   dims=3,
#                   perplexity=5,
#                   verbose=F,
#                   max_iter=500,
#                   check_duplicates=F)
# tsne <- data.frame(tSNE1 = tSNE_res[["Y"]][,1],
#                    tSNE2 = tSNE_res[["Y"]][,2],
#                    tSNE3 = tSNE_res[["Y"]][,3],
#                    cluster = cla)
# 
# tsne$colors <- ifelse(tsne$cluster %in% "GSE75010","#1D7E8F",ifelse(tsne$cluster %in% "GSE30186","#F0A274","#B8DADC"))
# colnames(tsne) <- c("x","y","z","cluster","colors")
# tsne$cluster <- as.factor(tsne$cluster)
# # pdf("4-PCA_After.pdf",width = 8.1,height = 6)
# pca_after <- scatterplot3d(tsne[,1:3],color = tsne$colors,main="After Normalization",pch = 21,bg = tsne$colors)
# legend("bottom",col = "black", legend = levels(tsne$cluster),pt.bg = c("#99cc00","#fcd300","#800080"), pch = 21,
#        inset = -0.2, xpd = TRUE, horiz = TRUE)
# dev.off()

#制作合并数据集的分组信息
adjusted_pd <- combined_pd
#对合并数据集的分组进行排序，Control在前，PE样本在后
adjusted_pd <- adjusted_pd[order(adjusted_pd$Group,decreasing = F),]
#保存合并数据集的分组信息
write.csv(adjusted_pd,"Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_pd,"../../3-DiffAnalysis/Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_pd,"../../6-WGCNA/Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_pd,"../../7&8-LASSO&SVM/Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_pd,"../../11-PE_ExpDiff&Cor//Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_pd,"../../12-PE_Immune/Combined_Datasets_Group.csv",row.names = F)

#根据合并数据集的分组信息对合并数据集样本进行重新排列
adjusted_gse<-adjusted_gse[,adjusted_pd$geo_accession] #对合并数据集样本进行重新排列
write.csv(adjusted_gse,"Combined_Datasets_Matrix.csv")
write.csv(adjusted_gse,"../../3-DiffAnalysis/Combined_Datasets_Matrix.csv")
write.csv(adjusted_gse,"../../6-WGCNA/Combined_Datasets_Matrix.csv")
write.csv(adjusted_gse,"../../7&8-LASSO&SVM/Combined_Datasets_Matrix.csv")
write.csv(adjusted_gse,"../../11-PE_ExpDiff&Cor/Combined_Datasets_Matrix.csv")
write.csv(adjusted_gse,"../../12-PE_Immune/Combined_Datasets_Matrix.csv")

write.csv(adjusted_gse[,adjusted_pd$Group == "PE"],"Disease_Matrix.csv")
write.csv(adjusted_gse[,adjusted_pd$Group == "PE"],"../../7&8-LASSO&SVM/Disease_Matrix.csv")
write.csv(adjusted_gse[,adjusted_pd$Group == "PE"],"../../9-Risk_Immune/Disease_Matrix.csv")
write.csv(adjusted_gse[,adjusted_pd$Group == "PE"],"../../10-Risk_Immune/Disease_Matrix.csv")
write.csv(adjusted_gse[,adjusted_pd$Group == "PE"],"../../13-ConsensusCluster/Disease_Matrix.csv")

table(adjusted_pd$Group)
# Control      PE 
# 91      94 
###END#####


###进行小鼠和人的基因转换######
#安装homologene这个R包
#install.packages('homologene')
# #加载homologene这个R包
# library(homologene)
# 
# adjustedgse <- as.data.frame(adjustedgse)
# genelist<-rownames(adjustedgse)
# #使用homologene函进行转换
# #@genelist是要转换的基因列表
# #@inTax是输入的基因列表所属的物种号，10090是小鼠
# #@outTax是要转换成的物种号，9606是人
# genelistid <- homologene(genelist, inTax = 10090, outTax = 9606)
# genelistid <- genelistid[,c(1,2)]
# colnames(genelistid) <- c("ID","Homo")
# genelistid <- as.data.frame(genelistid)
# gsehomo <- adjustedgse
# gsehomo$ID <- row.names(gsehomo)
# 
# gsehomo <- merge(genelistid,gsehomo, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
# gsehomo <- na.omit(gsehomo)
# gsehomo <- aggregate(x=gsehomo,by = gsehomo$Homo %>% list(),FUN = mean)#多个探针对应一个基因取平均
# row.names(gsehomo) <- gsehomo$Group.1
# gsehomo <- gsehomo[,-c(1,2,3)]
# 
# write.csv(gsehomo,"Combined_Datasets_Matrix_Homo.csv")
# 
# #NRGs
# nrgs <- data.table::fread("NRGs.csv",data.table = F)
# nrgslist <- nrgs$`Gene Symbol`
# #使用homologene函进行转换
# #@genelist是要转换的基因列表
# #10090是小鼠,9606是人
# nrgslistid <- homologene(nrgslist, inTax = 9606, outTax = 10090)
# nrgslistid <- nrgslistid[,c(1,2)]
# colnames(nrgslistid) <- c("Homo","Mus")
# 
# write.csv(nrgslistid,"NRGs_Mus.csv")

# prdegs <- readxl::read_xlsx("prdegs.xlsx")
prdegs <- data.table::fread("PRDEGs.csv",data.table = F)
dat_boxplot <- adjusted_gse[prdegs$x,]
dat_boxplot <- t(dat_boxplot)
dat_boxplot <- cbind(adjusted_pd$Group,dat_boxplot)
write.csv(dat_boxplot,"1-Boxplot.csv",row.names = F)














