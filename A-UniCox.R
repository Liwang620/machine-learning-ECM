###Uni&Multi Cox######
rm(list = ls())
# for(i in libr){
#   if (!requireNamespace(i, quietly = TRUE))
#     install.packages(i)
# }

for(i in c("survival", "forestplot", "magrittr")){
  library(i, character.only = T)
}

cox_extr <- function(fit){
  fit_summary <- summary(fit)
  dat_res <- data.frame(
    row.names = rownames(fit_summary$coef),
    p.value = signif(fit_summary$coef[,"Pr(>|z|)"]),
    mean = signif(fit_summary$coef[,"exp(coef)"]),
    lower = signif(fit_summary$conf.int[,"lower .95"]),
    upper = signif(fit_summary$conf.int[,"upper .95"]),
    coef = signif(fit_summary$coef[,"coef"]),
    check.rows = F
  )
  dat_res <- signif(dat_res,digits = 3)
  return(dat_res)
}


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
cox_analysis <- function(time, event, covariates, data, type=c("unicox","mulcox"), pvalue=0.05, save=TRUE, path="./", name="unicox.csv"){
  data <- na.omit(data)
  type <- match.arg(type)
  if(!all(data[, event] %in% c(0,1))) stop("event must be either 0 or 1")
  if(type == "unicox") {
    #对covariates中的每个因素都构建一个公式，默认使用OS.time和OS
    formulas <- sapply(covariates,
                       function(x) as.formula(paste0('Surv(', time, ",", event, ')~',  paste0("`", x, "`", sep = ""))))
    dat_res <- lapply(formulas, function(x){coxph(x, data = data)}) %>% #进行单因素cox回归
      lapply(cox_extr) %>% #提取单因素cox回归结果
      c(use.names = F) %>% 
      do.call(what = rbind) %>% #讲结果合并为数据框 
      dplyr::filter(p.value <= pvalue)
  } else {
    dat_res <- as.formula(paste('Surv(', time, ",", event, ')~',  paste(covariates, collapse = "+"))) %>% #构建公式
      coxph(data=data) %>% #进行多因素cox回归
      cox_extr #提取多因素cox回归结果
  }
  #保存结果
  if(save) {
    if(is.null(path)|is.null(name)){
      stop("No path or name")
    }else{
      write.csv(dat_res, paste0(path, "/", name))
    }
  }
  return(dat_res)
}



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
cox_plot <- function(dat_res, type=c("unicox", "mulcox"), boxcol="red", height=4, width=5, save=TRUE, path="./", name="unicox.pdf"){
  type <- match.arg(type)
  dat_res <- dat_res %>% 
    dplyr::mutate(GeneID = rownames(.), HR = paste(dat_res$mean, " [", dat_res$lower, " - ", dat_res$upper, "]", sep = "")) %>%
    dplyr::select("GeneID", "p.value", "HR", "mean", "lower", "upper")
  colnames(dat_res)[2] <- "P Value"
  labeltext <- rbind(colnames(dat_res)[1:3], dat_res[, 1:3])
  HR <- rbind(rep(NA, 3), dat_res[, 4:6])
  #画图
  p <- forestplot(labeltext = labeltext,             
                  HR,
                  zero = 1, 
                  lwd.ci = 1,#HR线的宽度
                    title=ifelse(type == "unicox", "Univariable Cox Regression Analysis",
                               "multivariable Cox regression analysis"),
                  colgap = unit(2, 'mm'),
                  xlab = "HR", 
                  txt_gp=fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8),
                                 title = gpar(cex = 1.2), cex = 1) ,#坐标轴和变量的字体大小
                  boxsize = 0.1, 
                  col=fpColors(box = boxcol, line = "black", zero = "black")#修改颜色
  ) 
  if(save) {
    if(is.null(path)|is.null(name)){
      stop("No path or name")
    }else{
      pdf(paste0(path, "/", name), height = height, width = width)
      print(p)
      dev.off()
    }
  }
  print(p)
}  


#sample
##读取数据

dat_tumor <- data.table::fread("Disease_Matrix.csv", data.table = F)
rownames(dat_tumor) <- dat_tumor[,1]
dat_tumor <- dat_tumor[,-1]
pd <- data.table::fread("Combined_Datasets_Clinical.csv", data.table = F)
# pd <- pd[pd$status == "HNSC",]
# colnames(pd)[1] <- "RNAseq"
genes <- data.table::fread("MRDEGs.csv", data.table = F)
dat_tumor <- dat_tumor[,pd$geo_accession]
dat_tumor <- t(dat_tumor[genes$x,])
data <- cbind(pd[,2:3],dat_tumor)
# data <- na.omit(data)

covariates <- c(colnames(data)[3:ncol(data)])

##cox分析
dat_res <- cox_analysis(time = "OS.time", event = "OS", covariates = covariates, data = data, type = "unicox", pvalue = 0.10,
                        save = T, path = "./", name = "1-UnicoxGenes.csv")

##可视化
cox_plot(dat_res, type = "unicox", boxcol = "#EE7785", height=8, width = 8.1, save = T, path = "./", name = "1-Unicox_ForestPlot.pdf")
dev.off()


##UniCox
# r <- c()
# library(progress)
# library(survival)
# library(survminer)
# 
# for(i in 5:356){
#   cox1=coxph(Surv(OS.time.y,OS)~unicoxinput[,i],data=unicoxinput)
#   r=rbind(r,c(colnames(unicoxinput)[i],coef(summary(cox1))[1,c(1,2,5)]))
#   #pb$tick()
#   Sys.sleep(1 / 100)
# }
# r <- data.frame(r)
# colnames(r)=c('probe','coef','HR','Pvalue')
# r1=subset(r,Pvalue<0.01) 
# cox_res <- r1$probe
# write.csv(cox_res,"0-LASSO_Input.csv")
# 
# riskscore <- readxl::read_xlsx("2-LASSO_RiskScore.xlsx")
# riskscore <- riskscore[,c(1:3,13)]
# riskscore$risk_level = ifelse(riskscore$lasso.risk.score >= median(riskscore$lasso.risk.score),"HighRisk","LowRisk")
# riskscore <- riskscore[order(riskscore$risk_level),]
# write.csv(riskscore,"LASSO_Group.csv",row.names = F)
# 
# libr <- c("survival", "forestplot", "dplyr", "magrittr")
