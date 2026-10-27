# 追加 Panel I: 核心抗菌基因时间轴折线,填满 3x3 右下角
# 复用主脚本已存对象需重算,故此脚本独立重建 I 并重拼图
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(stringr); library(tidyr)
  library(ggplot2); library(RColorBrewer)
})
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")

q3 <- read_excel("GSE235874_QC_BQC_Q3Norm_matrixonly.xlsx")
seg <- read_excel("GSE235874_segmentSummary_readTypes.xlsx")
expr <- as.data.frame(q3); rownames(expr) <- make.unique(as.character(expr$TargetName))
expr$TargetName <- NULL; expr <- as.matrix(expr); storage.mode(expr) <- "numeric"
logexpr <- log2(expr + 1)
parts <- str_split_fixed(colnames(expr), " \\| ", 3)
meta <- data.frame(col=colnames(expr), scan=parts[,1], roi=as.integer(parts[,2]), segment=parts[,3])
seg2 <- seg %>% transmute(scan=`Scan name`, roi=as.integer(`ROI name`), segment=`Segment name`, condition=as.character(Sample))
meta <- meta %>% left_join(seg2, by=c("scan","roi","segment"))
tp_levels <- c("control","CIP_3hr","CIP_12hr","CIP_24hr","CIP_48hr","CIP_72hr","CIP_5d","CIP_7d","CIP_14d","CIP_30d")
tp_levels <- tp_levels[tp_levels %in% unique(meta$condition)]
meta$condition <- factor(meta$condition, levels=tp_levels)

# 核心抗菌肽/防御基因(F 热图里 12-24hr 爆发的那批)
key_genes <- c("Reg3a","Reg3b","Reg3g","Reg2","Lyz2","Dmbt1","Lcn2","Hamp2")
key_genes <- intersect(key_genes, rownames(logexpr))

df <- logexpr[key_genes,,drop=FALSE] %>% t() %>% as.data.frame() %>%
  rownames_to_column("col") %>% left_join(meta[,c("col","condition")], by="col") %>%
  pivot_longer(-c(col,condition), names_to="gene", values_to="expr") %>%
  group_by(condition, gene) %>% summarise(mean=mean(expr), se=sd(expr)/sqrt(n()), .groups="drop") %>%
  mutate(tp_idx=as.integer(factor(condition, levels=tp_levels)))

BASE <- 13; TITLE <- 14
bold_theme <- theme(text=element_text(face="bold"), axis.text=element_text(face="bold",color="black"),
  axis.title=element_text(face="bold"), legend.text=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.title=element_text(face="bold",size=TITLE))

pI <- ggplot(df, aes(tp_idx, mean, color=gene)) +
  geom_line(linewidth=1.2) + geom_point(size=2.2) +
  scale_x_continuous(breaks=seq_along(tp_levels), labels=tp_levels, expand=expansion(mult=0.02)) +
  scale_color_brewer(palette="Set1") +
  labs(x=NULL, y="log2 expression", title="I. Antimicrobial peptide genes over time") +
  guides(color=guide_legend(ncol=2)) +
  theme_bw(base_size=BASE) + bold_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=BASE-3,face="bold"),
        legend.position="bottom", legend.text=element_text(size=BASE-4,face="bold"),
        legend.key.height=unit(9,"pt"), legend.title=element_blank(),
        legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"))
saveRDS(pI, "results/panels/pI.rds")
cat("Panel I 基因:", paste(key_genes, collapse=", "), "\n")
