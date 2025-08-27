###############KEGG_Pathview#################################
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")
DEG <- data.table::fread("MRDEGs_logFC.csv",data.table = F)


#####1. 添加ENTREZID列##########
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
genelist <- bitr(DEG$V1, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("V1"="SYMBOL"))

gene_diff = DEG$ENTREZID
#save(gene_diff,file = "ENTREZID.rda")

#load("ENTREZID.rda")

#####2. KEGG分析##########
R.utils::setOption( "clusterProfiler.download.method",'wininet')
kk <- enrichKEGG(gene         =  gene_diff,
                 organism ="hsa", #物种人类 hsa 小鼠mmu
                 keyType = "kegg",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.1)
# write.table(kk,file="1-KEGG.txt",sep="\t",quote=F,row.names = F)
#save(kk, file = "KEGG.Rdata")

#####整理##########
deg <- data.table::fread("MRDEGs_logFC.csv",data.table = F)
deg_up <- deg[deg$group%in%"up",]
deg_down <- deg[deg$group%in%"down",]
deg_diff <- rbind(deg_up,deg_down)
pvalue_diff <- deg_diff[,"adj.P.Val"]
names(pvalue_diff)=rownames(deg_diff)

#4. 绘制通路图
library("pathview")
#keggxls=read.table("KEGG.txt",sep="\t",header=T)
#write.csv(keggxls,file="KEGGXLS.csv",row.names = T,quote = F)

##4.1 单个样本通路图可视化
#options(download.file.method = "auto")
#debug(download.kegg)
# download.file("https://rest.kegg.jp/get/hsa00450/kgml", "./hsa00450.xml", quiet = T, mode = "wb", method = "auto")
# match.arg(method, c("auto", "internal", "wininet", "libcurl", "wget", "curl", "lynx"))
# rvest::read_html("https://rest.kegg.jp/get/hsa00450/kgml")

hsa04510 <- pathview(gene.data  = -log10(pvalue_diff),
                     pathway.id = "hsa04510",
                     species    = "hsa",
                     out.suffix = "pathview", 
                     limit      = list(gene= 10, cpd= 1))


hsa05418 <- pathview(gene.data  = -log10(pvalue_diff),
                     pathway.id = "hsa05418",
                     species    = "hsa",
                     out.suffix = "pathview", 
                     limit      = list(gene= 10, cpd= 1))

hsa05323 <- pathview(gene.data  = -log10(pvalue_diff),
                     pathway.id = "hsa05323",
                     species    = "hsa",
                     out.suffix = "pathview", 
                     limit      = list(gene= 10, cpd= 1))

hsa05205 <- pathview(gene.data  = -log10(pvalue_diff),
                     pathway.id = "hsa05205",
                     species    = "hsa",
                     out.suffix = "pathview", 
                     limit      = list(gene= 10, cpd= 1))

hsa00330 <- pathview(gene.data  = -log10(pvalue_diff),
                     pathway.id = "hsa00330",
                     species    = "hsa",
                     out.suffix = "pathview", 
                     limit      = list(gene= 10, cpd= 1))
##4.2 所有通路图可视化
for(i in keggxls$ID){
  pv.out <- pathview(gene.data = -log10(pvalue_diff), 
                     pathway.id = i, 
                     species = "hsa", 
                     out.suffix = "pathview",
                     limit=list(gene=10, cpd=10),
                     bins = list(gene = 10, cpd= 10))
}

