suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(stringr)
  library(limma); library(clusterProfiler); library(org.Mm.eg.db)
})
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")
# 重建 core 和 deg12(与主脚本一致)
q3 <- read_excel("GSE235874_QC_BQC_Q3Norm_matrixonly.xlsx")
seg <- read_excel("GSE235874_segmentSummary_readTypes.xlsx")
expr <- as.data.frame(q3); rownames(expr) <- make.unique(as.character(expr$TargetName))
expr$TargetName <- NULL; expr <- as.matrix(expr); storage.mode(expr)<-"numeric"; logexpr<-log2(expr+1)
parts <- str_split_fixed(colnames(expr)," \\| ",3)
meta <- data.frame(col=colnames(expr),scan=parts[,1],roi=as.integer(parts[,2]),segment=parts[,3])
seg2 <- seg %>% transmute(scan=`Scan name`,roi=as.integer(`ROI name`),segment=`Segment name`,condition=as.character(Sample))
meta <- meta %>% left_join(seg2,by=c("scan","roi","segment"))
tp <- c("control","CIP_3hr","CIP_12hr","CIP_24hr","CIP_48hr","CIP_72hr","CIP_5d","CIP_7d","CIP_14d","CIP_30d")
tp <- tp[tp %in% unique(meta$condition)]; meta$condition <- factor(meta$condition,levels=tp); meta<-meta[order(meta$condition),]
logexpr <- logexpr[,meta$col]
design <- model.matrix(~0+condition,data=meta); colnames(design)<-levels(meta$condition)
fit <- lmFit(logexpr,design); cips<-setdiff(levels(meta$condition),"control")
fit2 <- eBayes(contrasts.fit(fit,makeContrasts(contrasts=paste0(cips,"-control"),levels=design)))
dg <- list()
for(cc in cips){ tt<-topTable(fit2,coef=paste0(cc,"-control"),number=Inf); 
  dg[[cc]]<-rownames(tt)[tt$adj.P.Val<0.05 & abs(tt$logFC)>1] }
vl <- dg[c("CIP_3hr","CIP_12hr","CIP_24hr")]; core <- Reduce(intersect,vl)

cat("=== core (韦恩交集) 基因数:", length(core), "===\n")
m1 <- suppressWarnings(bitr(core,"SYMBOL","ENTREZID",org.Mm.eg.db))
fail1 <- setdiff(core, m1$SYMBOL)
cat("core 映射失败", length(fail1), "个 (", round(100*length(fail1)/length(core),2), "%):", paste(fail1,collapse=", "), "\n\n")

cat("=== deg12 (12hr全部DEG) 基因数:", length(dg[["CIP_12hr"]]), "===\n")
m2 <- suppressWarnings(bitr(dg[["CIP_12hr"]],"SYMBOL","ENTREZID",org.Mm.eg.db))
fail2 <- setdiff(dg[["CIP_12hr"]], m2$SYMBOL)
cat("deg12 映射失败", length(fail2), "个 (", round(100*length(fail2)/length(dg[["CIP_12hr"]]),2), "%):", paste(fail2,collapse=", "), "\n\n")

# 复现性检验:同一 core 连跑两次,失败率是否一致
a <- suppressWarnings(bitr(core,"SYMBOL","ENTREZID",org.Mm.eg.db))
b <- suppressWarnings(bitr(core,"SYMBOL","ENTREZID",org.Mm.eg.db))
cat("=== 复现性:同一core两次映射行数", nrow(a), "vs", nrow(b), "->",
    ifelse(nrow(a)==nrow(b),"一致(稳定)","不一致(有漂移!)"), "===\n")
