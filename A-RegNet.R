#基于小鼠基因的转录因子预测
library(magrittr)

## genes----
gene <- data.table::fread("cor_input.csv",data.table = F) %>% colnames()

## chipbase数据库----
## chipbase-tf download----
library("AnnotationDbi")
#BiocManager::install("org.Mm.eg.db")
library('org.Mm.eg.db')

ENSEMBL <- mapIds(org.Mm.eg.db,keys=gene,column="ENSEMBL",keytype="SYMBOL",multiVals="first")

for(i in 1:length(gene)){
  # i = 1
  gene_name <- gene[i]
  ensembl_name <- ENSEMBL[i]
  file = paste("chipbase/",gene_name,".xlsx",sep = "")
  link = paste( "https://rnasysu.com/chipbase3/download.php?base_page=regulator&organism=mouse&assembly=mm39&ref_gene_id=", ensembl_name, "&gene_symbol=", gene_name, 
                "&protein=0&regulator_type=tf&upstream=1kb&downstream=1kb&up_down_flag=0&motif_status=Y&sample_flag=0&protein_flag=0&Ftype=xls", sep = "")
  print(link)
  print(file)
  download.file(link,file)
}
all_files <- list.files("chipbase", full.names = T)
list_data <- list()

for(i in 1:length(all_files)){
  # i = 1
  message(Sys.time(), "  ", i)
  gene_name <- sub(".*?/(.*?)\\.xlsx$", "\\1", all_files[i])
  list_data[[gene_name]] <- data.table::fread(all_files[i], data.table = F) %>%
    dplyr::select(tf = "Transcription factor",upstream = "Number of binding sites contain motif (upstream)"
                  ,downstream = "Number of binding sites contain motif (downstream)") %>%
    dplyr::mutate(mrna = gene_name) %>%
    dplyr::distinct() ## 去重，按所有的列
}
chipbase <- do.call(rbind, list_data)
motifsum <- apply(chipbase[,2:3],1,sum)
table(motifsum)
tf_mrna <- chipbase[which(motifsum > 3),c(4,1)]

write.csv(tf_mrna,"3-mRNA-tf_node.csv",row.names = F, quote = F)

dat_id_type <- tf_mrna %>% 
  tidyr::gather(key = "type", value = "id") %>% 
  dplyr::select(2,1) %>% 
  dplyr::distinct()
openxlsx::write.xlsx(dat_id_type, file = "4-mRNA-tf_attribute.xlsx")
table_tf <- tf_mrna
colnames(table_tf)[1] <- "mRNA"
colnames(table_tf)[2] <- "TF"
write.csv(table_tf,"TableS1 mRNA-TF.csv",row.names = F,quote = F)


## mirna----
miRNA_Tarbase <- data.table::fread("TarBase_v8_download.txt",data.table = F)
# miRNA_Mus <- miRNA_Tarbase[miRNA_Tarbase$species == "Mus musculus",]
miRNA_Hs <- miRNA_Tarbase[miRNA_Tarbase$species == "Mus musculus",]
# miRNA_Hs <- miRNA_Hs[miRNA_Hs$geneName %in% hubgenes$hubGenes,]
# miRNA_Hs <- miRNA_Hs[miRNA_Hs$direct_indirect == "DIRECT",]
# mRNA_miRNA <- miRNA_Hs[which(is.na(miRNA_Hs$condition)),]
# mRNA_miRNA <- mRNA_miRNA[which(!is.na(mRNA_miRNA$tissue)),]
mRNA_miRNA <- miRNA_Hs[miRNA_Hs$geneName %in% gene,]
# mRNA_miRNA <- mRNA_miRNA[which(is.na(mRNA_miRNA$cell_line)),]
# mRNA_miRNA <- mRNA_miRNA[mRNA_miRNA$category == "Cancer/Malignant",]


mRNA_miRNA <- mRNA_miRNA[,c(2,3)]
colnames(mRNA_miRNA)[1] <- "node1"
colnames(mRNA_miRNA)[2] <- "node2"
write.csv(mRNA_miRNA,"1-mRNA-miRNA_nodes.csv",row.names = F,quote = F)
miRNA_attr <- data.frame("node" = c(unique(mRNA_miRNA$node1),unique(mRNA_miRNA$node2)),
                         "attribute" = c(rep("mRNA",length(unique(mRNA_miRNA$node1))),
                                         rep("miRNA",length(unique(mRNA_miRNA$node2)))))
write.csv(miRNA_attr,"1-mRNA-miRNA_attribute.csv",row.names = F,quote = F)

table_miRNA <- mRNA_miRNA
colnames(table_miRNA)[1] <- "mRNA"
colnames(table_miRNA)[2] <- "miRNA"
write.csv(table_miRNA,"TableS2 mRNA-miRNA.csv",row.names = F,quote = F)
