HLA <- data.table::fread("GeneCards-SearchResults.csv", data.table = F)
HLA <- as.data.frame(HLA)
TCGA_disase_mat <- data.table::fread('Disease_Matrix_FPKM.csv')
TCGA_disase_mat <- as.data.frame(TCGA_disase_mat)
rownames(TCGA_disase_mat) <- TCGA_disase_mat$V1
TCGA_disase_mat <- TCGA_disase_mat[,-1]
group <- as.data.frame(data.table::fread('Lasso-Riskscore.csv'))
TCGA_disase_mat <- TCGA_disase_mat[,group$sample_id]

HLA <- TCGA_disase_mat[HLA$`Gene Symbol`,]
HLA <- na.omit(HLA)
HLA <- t(HLA)
HLA <- as.data.frame(HLA)
HLA$group <- group$`Risk Level`

write.csv(HLA,'HLA_分组比较图输入.csv',row.names = FALSE)

#####免疫检查点#####
checkpoint_gene <- as.data.frame(data.table::fread("免疫检查点.csv"))

checkpoint_gene_mat <- TCGA_disase_mat[checkpoint_gene$Gene,]
checkpoint_gene_mat <- na.omit(checkpoint_gene_mat)

checkpoint_gene_mat <- t(checkpoint_gene_mat)
checkpoint_gene_mat <- as.data.frame(checkpoint_gene_mat)
checkpoint_gene_mat$trt <- group$`Risk Level`
write.csv(checkpoint_gene_mat,'免疫检查点分组比较图输入2.csv')

####ICD#####
ICD_gene <- as.data.frame(data.table::fread("ICD_genes.csv"))
ICD_gene_mat <- TCGA_disase_mat[ICD_gene$Genes,]
ICD_gene_mat <- na.omit(ICD_gene_mat)
ICD_gene_mat <- as.data.frame(t(ICD_gene_mat))
ICD_gene_mat$trt <- group$`Risk Level`
write.csv(ICD_gene_mat,'ICD分组比较图输入.csv')

