rm(list = ls())

###INPUT######
dat <- data.table::fread("Combined_Datasets_Matrix.csv",data.table = F)
row.names(dat) <- dat[,1]
dat <- dat[,-1]

# gene <- readxl::read_xlsx("hubGenes.xlsx")
gene <- data.table::fread("ModuleGenes.csv",data.table = F)

###CIBERSORT######
table_cell <- read.table('LM22_input.txt',fill=T,header=T,sep='\t',row.names=1,check.names=F)
data <- dat[intersect(rownames(dat),rownames(table_cell)),]
table_cell <- table_cell[intersect(rownames(data),rownames(table_cell)),]
source('A-CIBERSORT_Source.R')
result <- CIBERSORT(table_cell, data, perm = 1000, QN = T) 
# write.table(as.data.frame(result)[,1:22],"1-CIBERSORT_Results.txt",row.names = T,sep = "\t")

cibersort <- as.data.frame(result)[,1:22]
group <- data.table::fread("Combined_Datasets_Group.csv", data.table = F, header = T)
cibersort <- cibersort[group$geo_accession,] 

###可视化######
#成分柱状图
ciber_filter <- cibersort[,colSums(cibersort) > 0]   #留下丰度和大于0的细胞
conNum=nrow(group[group$Group=="Control",]) ## Control组数量
treatNum=nrow(group[group$Group==	"PI",]) ## Disease组数量
colA <- "#2486b9" # Control组颜色
colB <- "#f0a274" # Disease组颜色

pdf("1-Barplot.pdf",16.6,5)
set.seed(500) ## 种子
randomcoloR::distinctColorPalette(22) 
library(ggplot2)
mycol <- alpha(rainbow(ncol(ciber_filter)), 0.4) #创建彩虹色板（带60%透明度）
par(bty="o", mgp = c(2.5,0.3,0), ## bty边框类型，o就好；mgp：标题，坐标轴标题，坐标轴离边框距离
    mar = c(4.1,4.1,2.1,10.1), ## 空白边界行数(下左上右)
    tcl=-.25,las = 1,xpd = F) ## tcl坐标轴刻度线长度；las = 1：水平方向；xpd = F：打印剪裁到打印区域
a1 <- barplot(as.matrix(t(ciber_filter)),space = 0,
              names.arg = rep("",nrow(ciber_filter)), # 无横坐标样本名
              yaxt = "n", # 先不绘制y轴
              main = bquote(''),
              ylab = "Relative percentage", # 修改y轴名称
              col = mycol,
              border = NA)  # 修改是否包含外框
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # side = 2,坐标轴处于(下左上右)的左；补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2], # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber_filter), 
       xpd = T,fill = mycol,cex = 0.8, 
       border = NA, y.intersp = 1,
       x.intersp = 0.2,bty = "n")
par(srt=0,xpd=T)
rect(xleft = 0, ybottom = -0.01, xright = a1[conNum]+(a1[1]-0), ytop = -0.06,col=colA,border = "black")
text(a1[conNum]/2,-0.033,"Control",cex=1) ## Control组标签位置-0.033
rect(xleft = a1[conNum]+(a1[1]-0), ybottom = -0.01, xright =a1[length(a1)]+(a1[1]-0) , ytop = -0.06,col=colB,border = "black")
text((a1[length(a1)]+a1[conNum])/2,-0.033,"PE",cex=1) ## Disease组标签位置-0.033
dev.off() # 关闭画板

#######################################################
#绘制免疫细胞在不同分组中的分组比较图
#计算和筛选
boxplot_p <- function(data,vs,group){#vs：要进行比较的列名。
  formulas <- sapply(vs,
                     function(x) as.formula(paste(x,"~",group)))
  group <- data[[group]]
  if(length(levels(factor(group)))>2){
    test=lapply(formulas, function(x){kruskal.test(x, data = data)})
  }else{
    test=lapply(formulas, function(x){wilcox.test(x, data = data)})
  }
  test_results <- lapply(test,
                         function(x){ 
                           p.value<-signif(x$p.value, digits=3)
                           pstar <-ifelse(x$p.value<0.05,ifelse(x$p.value<0.01,"**","*"),"ns")
                           res <- c(p.value,pstar)
                           return(res)
                         })
  
  res <- t(as.data.frame(test_results, check.names = FALSE))
  colnames(res) <- c("p.value","pstar")
  return(res)
}

dat_boxplot <- as.data.frame(cibersort)
dat_boxplot <- cbind(group$Group,dat_boxplot)
colnames(dat_boxplot) <- gsub(" ","_",colnames(dat_boxplot))
colnames(dat_boxplot)[1] <- "Group"
colnames(dat_boxplot)[10] <- "T_cells_regulatory_Tregs"
dat_boxplot_p <- boxplot_p(dat_boxplot,colnames(dat_boxplot)[-1],"Group")

dat_boxplot_p <- as.data.frame(dat_boxplot_p)
dat_boxplot_p <- na.omit(dat_boxplot_p)
length(which(dat_boxplot_p$pstar != "ns"))
# [1] 11
paste0(rownames(dat_boxplot_p)[dat_boxplot_p$pstar != "ns"],collapse = "，")
# [1] "B_cells_naive，B_cells_memory，Plasma_cells，T_cells_CD8，T_cells_regulatory_Tregs，
# NK_cells_activated，Macrophages_M0，Macrophages_M2，Mast_cells_resting，Eosinophils，Neutrophils"

#可视化
ciber_boxplot <- ciber_filter[,dat_boxplot_p$pstar != "ns"]
immu<- rep(colnames(ciber_boxplot),each=nrow(ciber_boxplot)) #组别变量
immu <- factor(immu) #组别因子化
a<-c(group$Group)
boxplot_group <- rep(a,ncol(ciber_boxplot)) #每个组别的两个属性
boxplot_group <- factor(boxplot_group,levels = c('Control','PE')) #属性因子化
value <- c() #随机赋值
for (j in 1:ncol(ciber_boxplot)) { value<-c(value,ciber_boxplot[,j])}
value<-as.numeric(value)
boxplot_dat <- data.frame(immu_cell=immu,group=group,value=value) #生成数据框

library(ggthemes)
library(ggpubr)
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

boxplot <- 
  ggplot(boxplot_dat,aes(x=immu_cell,y=value,fill=boxplot_group)) +
  geom_boxplot(width=0.7*0.47,size=0.3*0.47,outlier.color = NA) +
  scale_fill_manual(values = c("PE" = "#f0a274", "Control" = "#2486b9")) +#颜色顺序同前面group因子化顺序一致。
  mytheme +
  stat_compare_means(method	= "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif",size = 2)+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position = "top",
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(),
        legend.text = element_text()) + # 隐藏图注标签 
  xlab('') +
  ylab('Infiltration Abundance')
boxplot
ggsave(boxplot,file="2-DiffBoxplot.pdf",units = "cm",width = 16.6,height = 7.0)
# output1 <- output1[-c(5,7),]
# rownames(output1)

#绘制免疫细胞相关性热图
heatmap <- t(ciber_filter[,dat_boxplot_p$pstar != "ns"])
# row.names(output[order(row.names(output)),])
#保留在分组比较图中全部细胞
# cells <- row.names(heatmap)
#immu<-t(gsva_es1[which(rownames(gsva_es1)%in% cells),which(group1$Group)])
heatmap <- t(heatmap[,group$geo_accession])
# heatmap <- t(heatmap)
# heatmap <- heatmap[,colSums(heatmap) > 0]
heatmap_col <- colorRampPalette(colors = c("#83a78d","white","#eea2a4"),space="Lab")#定义颜色,相关性排序：低，中，高
library(corrplot)
#cor <- cor(immu)
#method参数圆"circle", 方块"square", 椭圆"ellipse", 数字"number", 阴影"shade", 颜色"color", 饼形图"pie"
pdf("3-HeatMap.pdf",width = 15,height = 13)
corrplot(corr =cor(heatmap,method = 'spearman'),order="AOE",type="upper",tl.pos="tp",method="color",
         tl.col = "black",col=heatmap_col(50))
corrplot(corr = cor(heatmap,method = 'spearman'),add=TRUE, type="lower", method="number",order="AOE", 
         diag=FALSE,tl.pos="n", cl.pos="n",col="black")#col1(50))
dev.off()
#绘制在High分组中的免疫细胞相关性热图。High样本87例
#immu<-t(gsva_es1[which(rownames(gsva_es1)%in% cells),which(group1$Group == "Control")])
#library(corrplot)
#pdf("4-High_HeatMap.pdf",width = 12,height = 10)
#corrplot(corr =cor(immu,method = 'spearman'),order="AOE",type="upper",tl.pos="tp",method="color",
#         tl.col = "black",col=col1(50))
#corrplot(corr = cor(immu,method = 'spearman'),add=TRUE, type="lower", method="number",order="AOE",
#         diag=FALSE,tl.pos="n", cl.pos="n",col="black")#col1(50))
#dev.off()
#绘制22个免疫细胞与8目的基因的相关性热图。
hub <- gene$x
corheatmap_hub <- dat[,group$geo_accession]
corheatmap_hub <- corheatmap_hub[hub,]

#合并免疫细胞样本矩阵与目的基因的样本矩阵，前21列为免疫细胞（1:8），后面6个为目的基因（9:11）
corheatmap_dat <- as.data.frame(ciber_filter[,dat_boxplot_p$pstar != "ns"])
# corheatmap_dat <- corheatmap_dat[,colSums(corheatmap_dat) > 0]
corheatmap <- cbind(corheatmap_dat,t(corheatmap_hub))

library("Rmisc")
library("plyr")
library("Hmisc")
cor<-data.frame()
for (mm in 12:15) {#目的基因所在列
  cor1<-data.frame(0,0,0,0)
  for (i in 1:11) {#免疫细胞所在列
    c<-rcorr(corheatmap[,i],corheatmap[,mm],type = 'spearman')
    cor1[i,1]<-c$r[2]
    cor1[i,2]<-c$P[2]
    cor1[i,3]<-colnames(corheatmap)[i]
    cor1[i,4]<-colnames(corheatmap)[mm]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("correlation","Pvalue","immu_cell","gene")
cor_filter <- cor[cor$Pvalue<0.05,]
cor_filter<-na.omit(cor_filter)
cor_filter<-cor_filter[order(cor_filter$correlation,decreasing = F),]
y=factor(cor_filter$gene)
x=factor(cor_filter$immu_cell)
a<-(-log10(cor_filter$Pvalue))

# library(corrplot)
corheatmap <- 
  ggplot(cor_filter,aes(x,y)) + 
  geom_point(aes(size=a,color=correlation)) +
  scale_colour_gradient2(high="#eea2a4",mid = "white",low="#83a78d") + 
  mytheme +
  # theme(plot.title = element_text(hjust = 0.5)) +
  labs(size="-log10Pvalue",x="",y="") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key = element_blank(),
        legend.margin = margin(l = 0.2,unit = "cm"))
corheatmap
ggsave(corheatmap,file="4-CorHeatMap.pdf",unit="cm", width =12,height =10)
###END######


