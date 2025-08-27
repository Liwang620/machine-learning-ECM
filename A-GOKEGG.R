rm(list=ls())   #清空环境
options(stringsAsFactors = F)#设置变量非因子化

# 数据处理 --------------------------------------------------------------------

DEG <- data.table::fread("input-DEGs.txt",data.table = F)#读入数据源
#只需要基因名称，logFC两列信息。
colnames(DEG)#查看DEG表格的列名
DEG <- DEG[,c("ID","logFC")] #仅留下用到的列
## #加载对应的包，先添加ENTREZID列
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
#获得基因对应的ENTREZID信息
genelist <- bitr(DEG$ID, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("ID"="SYMBOL"))#在表格中加入ENTREZID列

DEG$group <-ifelse(DEG$logFC >= 1,"Up","Down")#根据logFC的值添加上下调的信息。


gene_up = DEG[DEG$group == 'Up','ENTREZID'] # 上调基因
gene_down = DEG[DEG$group == 'Down','ENTREZID'] # 下调基因
gene_diff = c(gene_up,gene_down) # 把全部基因合并起来
write.table(DEG,file="1-DEGs-new.txt",sep="\t",quote=F,row.names = F)#保存表格


#1.GO富集分析
GO1<-enrichGO(gene_diff,#输入基因对应的ENTREZID
              OrgDb ='org.Hs.eg.db',#参数 OrgDb 指定该物种对应的 org 包的名字，
              pAdjustMethod="BH",#pAdjustMethod 指定多重假设检验矫正的方法，
              ont = "ALL",#ont 代表 GO 的 3 大类别，BP, CC, MF，也可以是全部的 all
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.1,#设定q值阈值
              readable = T)#readable=TRUE代表将基因 ID 转换为 gene symbol。

GO1@result=GO1@result[order(GO1@result$Count,decreasing = T),]#按通路中包含的输入基因数量从高到底排序

Go_res <- GO1@result#把结果转成表格
write.csv(Go_res,file="2-GO-12DEGs.csv",row.names = T,quote = F)#保存GO分析的结果
BP <- subset(Go_res,ONTOLOGY %in% "BP")#挑选并保存结果中的BP通路
CC <- subset(Go_res,ONTOLOGY %in% "CC")#挑选并保存结果中的CC通路
MF <- subset(Go_res,ONTOLOGY %in% "MF")#挑选并保存结果中的MF通路

Go_outcome<- rbind(BP[1:5,],CC[1:5,],MF[1:5,])#GO通路一般挑选不出来明星通路，所以按基因数量选三种通路的前五个，也可以自己在外部筛选。
Go_outcome <- na.omit(Go_outcome)#去除空行
rownames(Go_outcome)#查看结果表格的行名


#####2. KEGG分析##########
R.utils::setOption( "clusterProfiler.download.method",'wininet')#设置clusterProfiler包内部的一些标准
kk <- enrichKEGG(gene         =  gene_diff,#输入基因对应的ENTREZID
                 organism ="hsa", #物种人类 hsa 小鼠mmu
                 keyType = "kegg",#分析类型"kegg"
                 pAdjustMethod = "BH",#调整方法选择BH
                 pvalueCutoff = 0.05,#阈值筛选标准，P值
                 qvalueCutoff =0.1)#阈值筛选标准，q值
                 
kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")#把kk结果中的ENTREZID转换成基因
kk1 <- kk@result#如果kk1没有表格，说明设置的标准里没富集到通路，不用往下做了。
#保存kegg分析结果
write.table(kk1,file="3-KEGG-12DEGs.txt",sep="\t",quote=F,row.names = F)#保存kegg分析结果

write.csv(kk1,file="4-KEGG-12DEGs.csv",row.names = T,quote = F)
#外部打开文件，"4-KEGG-12DEGs.csv"，仅保留想要展示的通路，然后保存。

#读取保存的通路
KEGGDEGs <- read.csv("4-KEGG-12DEGs.csv",header = T)#读取保存下来的KEGG文件

kk1 <- kk1[KEGGDEGs$ID,]#仅保留画图的通路
kk1$ONTOLOGY <- "KEGG"#添加ONTOLOGY列
#合并GO,KEGG的富集分析结果
Go_Kegg <- rbind(Go_outcome,kk1)

rownames(Go_Kegg)#查看结果表格的行名
Golevels <- as.character(Go_Kegg$ONTOLOGY[which(!duplicated(Go_Kegg$ONTOLOGY))])#查看分析结果的类别

Go_Kegg$ONTOLOGY <- factor(Go_Kegg$ONTOLOGY,levels= Golevels)#把分析的类别因子化，对分析的类别进行分类

Go_Kegg$type_order <- factor(Go_Kegg$Description,levels= Go_Kegg$Description)#通路因子化
# GO可视化
#画富集结果的柱状图
library(ggplot2)
ggplot(data = Go_Kegg)+
  geom_bar(aes(x = type_order,y = -log10(p.adjust),fill = ONTOLOGY),stat = "identity")+facet_wrap(~ ONTOLOGY, scales = "free")+#根据ONTOLOGY进行分面
  scale_fill_manual(values = c("#576584","#96222D","#00A087","#3C5488"))+#颜色顺序是类别的顺序
  #coord_flip()+#翻转坐标轴
  xlab("Term")+
  theme_classic()+
  theme(axis.text.x=element_text(face = "bold",angle = 40,vjust = 1, hjust = 1 ))
ggsave("5-GO_KEGG.pdf",width = 16,height = 7)#保存富集柱状图结果

#画富集结果的气泡图

ggplot(data = Go_Kegg,aes(x=type_order,y = Count,colour=pvalue))+
  geom_point(aes(size=Count),shape=16,stat = "identity")+
  scale_color_gradient(low = "#576584",high = "#96222D")+#设置颜色
  facet_wrap( ~ ONTOLOGY, scales = "free")+#根据ONTOLOGY进行分面
  scale_fill_manual(values = c("#576584","#96222D","#00A087","#3C5488"))+#颜色顺序是类别的顺序
  #coord_flip()+#翻转坐标轴
  xlab("Term")+
  theme_classic()+
  theme(axis.text.x=element_text(face = "bold",angle = 40,vjust = 1, hjust = 1 ))
ggsave("6-GO_KEGG-point.pdf",width = 16,height = 7)#保存富集柱状图结果


#制作logFC表格
library(enrichplot)
List = DEG$logFC
names(List)= DEG$ENTREZID
List = sort(List,decreasing = T)
# GO-BP网络图
showBP <- Go_Kegg$Description[which(Go_Kegg$ONTOLOGY == "BP")]#挑选BP通路
pdf(file="7-GOcircos-BP.pdf", width=11, height=7)
cnetplot(GO1, #数据使用GO分析得到的最初表格
                  categorySize="pvalue", 
                  showCategory = showBP,#只展示挑选的BP通路
                  foldChange = List,#使用整理好的logFC值
                  circular = TRUE,#画成环状，不想画成环状，直接注释这一行
                  node_label = "all",#展示图中全部信息
                  colorEdge = TRUE)
dev.off()

#这次没有富集到CC通路，所以这里不画CC通路的网络图了
# GO-MF网络图
showMF <- Go_Kegg$Description[which(Go_Kegg$ONTOLOGY == "MF")]#挑选MF通路
pdf(file="8-GOcircos-MF.pdf", width=9, height=7)
cnetplot(GO1, #数据使用GO分析得到的最初表格
         categorySize="pvalue", 
         showCategory = showMF,#只展示挑选的MF通路
         foldChange = List,#使用整理好的logFC值
         circular = TRUE,#画成环状，不想画成环状，直接注释这一行
         node_label = "all",#展示图中全部信息
         colorEdge = TRUE)
dev.off()

# KEGG网络图
showKEGG <- Go_Kegg$Description[which(Go_Kegg$ONTOLOGY == "KEGG")]#挑选KEGG通路
#下面画图如果报错，就调整阈值重新进行KEGG分析
library(clusterProfiler)
kk <- enrichKEGG(gene         =  gene_diff,#输入基因对应的ENTREZID
                 organism ="hsa", #物种人类 hsa 小鼠mmu
                 keyType = "kegg",#分析类型""
                 pAdjustMethod = "BH",#调整方法选择BH
                 pvalueCutoff = 1,#阈值筛选标准，P值
                 qvalueCutoff =1)#阈值筛选标准，q值
                

kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")#把kk结果中的ENTREZID转换成基因
library(enrichplot)
pdf(file="9-GOcircos-KEGG.pdf", width=9, height=7)

cnetplot(kk, #数据使用KEGG分析得到的最初表格
         categorySize="pvalue", 
         showCategory = showKEGG,#只展示挑选的KEGG通路
         foldChange = List,#使用整理好的logFC值
         circular = TRUE,#画成环状，不想画成环状，直接注释这一行
         node_label = "all",#展示图中全部信息
         colorEdge = TRUE)

dev.off()






#制作画气泡图的表格
kegg=as.data.frame(kk)#转化成矩阵
kegg$ONTOLOGY <- "KEGG"#添加ONTOLOGY列
GO=as.data.frame(GO1)#转化成矩阵
Go_Kegg2 <- rbind(GO,kegg)#联合两个矩阵
#提取有用的列
go=data.frame(Category=Go_Kegg2$ONTOLOGY, ID=Go_Kegg2$ID, Term=Go_Kegg2$Description, Genes = gsub("/", ", ", Go_Kegg2$geneID), adj_pval = Go_Kegg2$p.adjust)
genelist <- data.frame(ID = DEG$ID, logFC = DEG$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)#制作出用于画图的表格
###画bubble图

# library(GOplot),不好用，放弃了
# pdf(file="10-GO-kegg bubble.pdf", width=8, height=6)
# p <- GOBubble(circ, title="Bubble plot",##设置标题
#          colour = c("#576584","#96222D","#00A087","#3C5488"),#颜色顺序是类别的顺序
#          labels = 3,##-log(adj p-value)>3
#          ID = TRUE,##TRUE显示符合标准的GO term ID，FALSE显示GO term name
#          table.legend = FALSE, table.col = FALSE, ##右侧表格设置，这一行改为TRUE时是加通路表格，颜色
#          bg.col = FALSE,#背景颜色设置
#          display = 'multiple'
# )
# p + scale_size_continuous(range = c(1,5)) +
#   theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
# dev.off()
# 




###画bubble图
#加载包
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
pdf(file="10-GObubble.pdf", width=8, height=6)
circ1 <- circ[!duplicated(circ$ID),]#去除重复的通路名
circ1levels <- as.character(circ1$category[which(!duplicated(circ1$category))])#查看分析结果的类别

circ1$category <- factor(circ1$category,levels= circ1levels)#把分析的类别因子化，对分析的类别进行分类
ggplot(circ1,aes(x=zscore,y=-log10(adj_pval)))+
  geom_point(aes(size=count,color=category),alpha=0.6)+
  scale_color_brewer(palette = "Accent")+
  theme_bw()+
  theme(
    legend.position = c("none")
  )+
  geom_text_repel(
    data = circ1[-log10(circ1$adj_pval)>-2,],
    aes(label = circ1[-log10(circ1$adj_pval)>-2,]$ID),#添加符合标准的ID
    size = 3, max.overlaps = 20,
    segment.color = "black", show.legend = FALSE )+
  facet_grid(.~category)#根据通路类别进行分面
dev.off()

#结束

