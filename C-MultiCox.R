###Uni&Multi Cox######
# dat_multicox <- dat_tumor[,prognosis_lasso$genes]
# dat_multicox <- dat_multicox[pd$geo_accession,]
# dat_multicox <- cbind(pd[,4:6],dat_multicox)

# grep("Not Available", dat_multicox$stage_event.tnm_categories.clinical_categories.clinical_M)
# grep("Not Available", dat_multicox$stage_event.tnm_categories.clinical_categories.clinical_N)
# grep("Not Available", dat_multicox$stage_event.tnm_categories.clinical_categories.clinical_T)
# grep("Discrepancy", dat_multicox$stage_event.tnm_categories.clinical_categories.clinical_T)
# dat_multicox <- dat_multicox[-c(75,144,383,384,385),]
# write.csv(dat_multicox,"../7-PrognosticAnalysis/1-MultiCox_Input.csv",row.names = F)


# dat_km <- data.table::fread("RiskScore.csv",data.table = F)
# pd <- pd[pd$RNAseq %in% dat_km$V1,]
# dat_km <- cbind(pd[,4:5],dat_km[,-1])
# write.csv(dat_km[,-4],"KM.csv")
# rm(list = ls())
# for(i in libr){
#   if (!requireNamespace(i, quietly = TRUE))
#     install.packages(i)
# }

# for(i in c("survival", "forestplot", "magrittr")){
#   library(i, character.only = T)
# }

# cox_extr <- function(fit){
#   fit_summary <- summary(fit)
#   dat_res <- data.frame(
#     row.names = rownames(fit_summary$coef),
#     p.value = signif(fit_summary$coef[,"Pr(>|z|)"]),
#     mean = signif(fit_summary$coef[,"exp(coef)"]),
#     lower = signif(fit_summary$conf.int[,"lower .95"]),
#     upper = signif(fit_summary$conf.int[,"upper .95"]),
#     coef = signif(fit_summary$coef[,"coef"]),
#     check.rows = F
#   )
#   dat_res <- signif(dat_res,digits = 3)
#   return(dat_res)
# }


#' 进行cox分析
#' 
#' @param time(numeric) 生存时间
#' @param event(factor) 生存结局, Alive=0, Dead =1，需要转换为0, 1
#' @param covariates(character) 自变量
#'   \cr 如果有age，age需要转变为numeric
#'   \cr 如果有stage，stage需要将IA, IB, IIA, IIB, ..之类的修改为I, II, III..
#' @param data(data.frame) 输入数据，包含cox所需的变量
#' @param type(character) 单因素cox选择"unicox"，多因素cox选择"mulcox"
#' @param save(logic) TRUE or FALSE, 是否需要保存
#' @param path(character) 保存路径, 默认当前文件夹
#' @param name(character) 文件命名
#' 
#' @return (data.frame) cox分析结果
#' 
#' @export
# cox_analysis <- function(time, event, covariates, data, type=c("unicox","mulcox"), pvalue=0.05, save=TRUE, path="./", name="unicox.csv"){
#   data <- na.omit(data)
#   type <- match.arg(type)
#   if(!all(data[, event] %in% c(0,1))) stop("event must be either 0 or 1")
#   if(type == "unicox") {
#     #对covariates中的每个因素都构建一个公式，默认使用OS.time和OS
#     formulas <- sapply(covariates,
#                        function(x) as.formula(paste0('Surv(', time, ",", event, ')~',  paste0("`", x, "`", sep = ""))))
#     dat_res <- lapply(formulas, function(x){coxph(x, data = data)}) %>% #进行单因素cox回归
#       lapply(cox_extr) %>% #提取单因素cox回归结果
#       c(use.names = F) %>%
#       do.call(what = rbind) %>% #讲结果合并为数据框
#       dplyr::filter(p.value <= pvalue)
#   } else {
#     dat_res <- as.formula(paste('Surv(', time, ",", event, ')~',  paste(covariates, collapse = "+"))) %>% #构建公式
#       coxph(data=data) %>% #进行多因素cox回归
#       cox_extr #提取多因素cox回归结果
#   }
#   #保存结果
#   if(save) {
#     if(is.null(path)|is.null(name)){
#       stop("No path or name")
#     }else{
#       write.csv(dat_res, paste0(path, "/", name))
#     }
#   }
#   return(dat_res)
# }



#' cox分析可视化 
#' @param data_res(data.frame) cox_analysis函数得到的数据框
#' @param type(character) 单因素cox选择"unicox"，多因素cox选择"mulcox"
#' @param boxcol(character) HR平均值方块的颜色
#' @param height, width(numeric) 图片的高度, 图片的宽度;尺寸按照pdf()标准
#' @param save(logic) TRUE or FALSE, 是否需要保存
#' @param path(character) 保存路径，默认当前文件夹
#' @param name(character) 图片的名字
#' 
#' @return 森林图
#' 
#' @export
# cox_plot <- function(dat_res, type=c("unicox", "mulcox"), boxcol="red", height=4, width=5, save=TRUE, path="./", name="unicox.pdf"){
#   type <- match.arg(type)
#   dat_res <- dat_res %>% 
#     dplyr::mutate(GeneID = rownames(.), HR = paste(dat_res$mean, " [", dat_res$lower, " - ", dat_res$upper, "]", sep = "")) %>%
#     dplyr::select("GeneID", "p.value", "HR", "mean", "lower", "upper")
#   colnames(dat_res)[2] <- "P Value"
#   labeltext <- rbind(colnames(dat_res)[1:3], dat_res[, 1:3])
#   HR <- rbind(rep(NA, 3), dat_res[, 4:6])
#   #画图
#   p <- forestplot(labeltext = labeltext,             
#                   HR,
#                   zero = 1, 
#                   lwd.ci = 1,#HR线的宽度
#                     title=ifelse(type == "unicox", "Univariable Cox Regression Analysis",
#                                "multivariable Cox regression analysis"),
#                   colgap = unit(2, 'mm'),
#                   xlab = "HR", 
#                   txt_gp=fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8),
#                                  title = gpar(cex = 1.2), cex = 1) ,#坐标轴和变量的字体大小
#                   boxsize = 0.1, 
#                   col=fpColors(box = boxcol, line = "black", zero = "black")#修改颜色
#   ) 
#   if(save) {
#     if(is.null(path)|is.null(name)){
#       stop("No path or name")
#     }else{
#       pdf(paste0(path, "/", name), height = height, width = width)
#       print(p)
#       dev.off()
#     }
#   }
#   print(p)
# }  


#sample
##读取数据

# dat_tumor <- data.table::fread("Matrix_Disease_FPKM.txt", data.table = F)
# rownames(dat_tumor) <- dat_tumor[,1]
# dat_tumor <- dat_tumor[,-1]
# pd <- data.table::fread("PD.csv", data.table = F)
# pd <- pd[pd$status == "Tumor",]
# mdrdegs <- data.table::fread("MDRDEGs.csv", data.table = F)
# dat_tumor <- dat_tumor[,pd$RNAseq样本编号]
# dat_tumor <- t(dat_tumor[mdrdegs$x,])
# data <- cbind(pd[,7:8],dat_tumor)
# data <- na.omit(data)

covariates <- c(prognosis_lasso$genes)

##cox分析
#注意：多因素cox事件数（OS为0的数目和为1的数目）一般要大于10倍的变量数（基因数）
dat_res_multi <- cox_analysis(time = "OS.time", event = "OS", covariates = covariates, data = data, type = "mulcox", pvalue = 0.05,
                        save = T, path = "./", name = "3-MulticoxGenes.csv")

##可视化
cox_plot(dat_res_multi, type = "mulcox", boxcol = "#EE7785", height=8, width = 8.1, save = T, path = "./", name = "4-Multicox_ForestPlot.pdf")
dev.off()

#RiskScore
#' riskscore计算
#' 
#' @param dat_res(data.frame) lasso_analysis分析结果
#' @param data(data.frame) 包含筛选基因的数据框
#' 
#' @return (data.frame) 包含riskscore和risklevel的数据框
#' 
#' @export
risk_score_mulcox <- function(dat_res_multi,dat_riskscore){
  genes <- rownames(dat_res_multi)
  multicox_coef <- dat_res_multi[, "coef"]
  riskscore <- dat_riskscore[, genes] %>% apply(1, function(x) {crossprod(as.numeric(x), multicox_coef)})
  risklevel <- ifelse(riskscore>=median(riskscore), "High Risk", "Low Risk")
  risk <- data.frame(riskscore, risklevel)
  return(risk)
}

riskscore_mulcox <- risk_score_mulcox(dat_res_multi,dat_riskscore)
riskscore_mulcox <- riskscore_mulcox[order(riskscore_mulcox$risklevel),]
riskscore_mulcox$sampleid <- rownames(riskscore_mulcox)
write.csv(riskscore_mulcox,"RiskScore_Mulcox.csv")
saveRDS(dat_res_multi, file = "../10-ModelValidation/RiskScore/Mulcox_Coef.rds")
# riskgroup <- readxl::read_xlsx("KM_Group.xlsx")
# colnames(riskgroup)[3] <- "riskscore"
# riskscore_mulcox <- cbind(riskscore_mulcox[,c(3,1)],riskgroup[,c(1,2,4)])
# riskscore_mulcox <- riskscore_mulcox[order(riskscore_mulcox$group),]
# riskscore_mulcox <- riskscore_mulcox[1:501,]

#RISK PLOT
dat_riskplot_mulcox <- dat_riskscore[riskscore_mulcox$sampleid,]
dat_riskplot_mulcox <- cbind(pd[pd$geo_accession %in% rownames(dat_riskplot_mulcox),2:3],riskscore_mulcox[,1:2],dat_riskplot_mulcox[,prognosis_lasso$genes])
rownames(dat_riskplot_mulcox) <- riskscore_mulcox$sampleid
# write.csv(dat_riskplot_mulcox,"../10-Risk_DiffAnalysis/1-RiskPlot.csv",row.names = T)
# write.csv(dat_riskplot_mulcox[,4,drop = F],"../10-Risk_DiffAnalysis/RiskGroup.csv",row.names = T)

dat_mulcox_score <- cbind(pd[pd$geo_accession %in% rownames(dat_riskplot_mulcox),-1],riskscore_mulcox[,1])
dat_mulcox_level <- cbind(pd[pd$geo_accession %in% rownames(dat_riskplot_mulcox),-1],riskscore_mulcox[,2])
write.csv(dat_mulcox_score,"../9-PrognosticAnalysis/0-Cox_Input_MulcoxScore.csv",row.names = F)
write.csv(dat_mulcox_level,"../9-PrognosticAnalysis/0-Cox_Input_MulcoxLevel.csv",row.names = F)

write.csv(dat_riskplot_mulcox,"../11-Risk_DiffAnalysis//1-Boxplot.csv",row.names = F)
# write.csv(dat_riskplot_mulcox[,-(1:4)],"../11-Risk_ExpDiff&ROC/3-Corheatmap.csv",row.names = F)
# dat_boxplot <- cbind(dat_riskplot_mulcox[,5],dat_riskplot_mulcox[,1:4])
# write.csv(dat_boxplot[,-1],"3-LIHC_Cor.csv",row.names = F)
#画LASSO模型的Nomogram
# library(rms)
# dc <- datadist(dat_riskplot_mulcox)
# options(datadist="dc")
# fit <- lrm(risklevel~riskscore, data=dat_riskplot_mulcox, x=T, y=T,tol=1e-9,maxit=1000)
# nom <- nomogram(fit,fun=plogis,fun.at =c(0.5,1) ,funlabel=c("Risk of LIHC"))
# #pdf(width = 8.1,height = 5,file = "Nomogram_1.pdf",onefile = F)
# plot(nom,cex.axis=1.5,cex.var=1.5)
# dev.off()

#画模型基因的Nomogram
# colnames(dat_riskplot_mulcox)
# fit <- lrm(risklevel~BIRC5+PLK1+NQO1+UCHL1+HMGA2, 
#            data=dat_riskplot_mulcox, x=T, y=T,tol=1e-9,maxit=1000)
# dc<-datadist(dat_riskplot_mulcox)
# options(datadist="dc")
# nom <- nomogram(fit,fun=plogis,fun.at =c(0.5,1) ,funlabel=c("Risk of LIHC"))
# pdf(width = 8.1,height = 5,file = "4-Nomogram.pdf",onefile = F)
# plot(nom,cex.axis=1.5,cex.var=1.5)
# dev.off()
