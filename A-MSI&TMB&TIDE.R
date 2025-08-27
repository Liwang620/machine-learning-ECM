rm(list = ls())


#TIDE
tide <- data.table::fread("TIDE.csv", data.table = F)
group <- data.table::fread("LASSO_Group.csv", data.table = F)
intersect_tide <- intersect(tide$Patient,group$sample_id)
tide <- tide[tide$Patient %in% intersect_tide,]
group_tide <- group[group$sample_id %in% intersect_tide,]
tide <- tide[group_tide$sample_id,]
boxplot_tide <- cbind(group_tide$risk_level,tide$TIDE)
boxplot_tide <- as.data.frame(boxplot_tide)
colnames(boxplot_tide)[1] <- "Group"
colnames(boxplot_tide)[2] <- "TIDE"
write.csv(boxplot_tide,"1-TIDE.csv",row.names = F)

#msi_tmb&TMB
msi_tmb_tmb <- data.table::fread("msi_tmb&TMB.tsv", data.table = F)
msi_tmb <- msi_tmb[!duplicated(msi_tmb$`Patient ID`),]
msi_tmb <- as.data.frame(msi_tmb)
rownames(msi_tmb) <- msi_tmb$`Patient ID`
group_msitmb <- group
group_msitmb$`Patient ID` <- stringr::str_sub(group_msitmb$sample_id,1,12)
boxplot_msitmb <- merge(msi_tmb,group_msitmb, by = "Patient ID")
boxplot_msi <- boxplot_msitmb[,c(66,30)]
boxplot_tmb <- boxplot_msitmb[,c(66,57)]
boxplot_msi <- boxplot_msi[order(boxplot_msi$risk_level),]
boxplot_tmb <- boxplot_tmb[order(boxplot_tmb$risk_level),]
colnames(boxplot_msi)[1] <- "Group"
colnames(boxplot_tmb)[1] <- "Group"
colnames(boxplot_msi)[2] <- "MSI"
colnames(boxplot_tmb)[2] <- "TMB"
boxplot_msi$MSI <- log2(boxplot_msi$MSI+1)
boxplot_tmb$TMB <- log2(boxplot_tmb$TMB+1)
write.csv(boxplot_msi,"2-MSI.csv",row.names = F)
write.csv(boxplot_tmb,"3-TMB.csv",row.names = F)

