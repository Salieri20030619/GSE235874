# ============================================================
# GSE235874 GeoMx WTA | 小鼠胰腺 CIP 急性胰腺炎时间轴
# 炎症通路 + response to bacterium(GO)| 参考 UC 图三 panel
# ============================================================
suppressPackageStartupMessages({
  need_cran <- c("readxl","dplyr","tibble","stringr","ggplot2","VennDiagram","ggalluvial","openxlsx")
  need_bioc <- c("limma","clusterProfiler","org.Mm.eg.db","enrichplot")
  miss_c <- need_cran[!need_cran %in% rownames(installed.packages())]
  miss_b <- need_bioc[!need_bioc %in% rownames(installed.packages())]
  if(length(miss_c)) install.packages(miss_c, repos="https://cloud.r-project.org")
  if(length(miss_b)){ if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(miss_b, update=FALSE, ask=FALSE) }
  library(readxl); library(dplyr); library(tibble); library(stringr)
  library(limma); library(ggplot2); library(VennDiagram); library(ggalluvial)
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot); library(openxlsx)
})
flog.threshold(ERROR, name="VennDiagramLogger")  # 关掉 VennDiagram 冗余日志
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")
dir.create("results", showWarnings=FALSE)

# ---- 参数 ----
LFC <- 1.0        # |log2FC| 阈值
PADJ <- 0.05      # 校正 p 阈值
q3_file  <- "GSE235874_QC_BQC_Q3Norm_matrixonly.xlsx"
seg_file <- "GSE235874_segmentSummary_readTypes.xlsx"

# ============ 1. 读入 + 构建 ROI→condition 元数据 ============
q3  <- read_excel(q3_file)
seg <- read_excel(seg_file)

expr <- as.data.frame(q3)
rownames(expr) <- make.unique(as.character(expr$TargetName))
expr$TargetName <- NULL
expr <- as.matrix(expr)                 # genes x ROI (Q3 线性值)
storage.mode(expr) <- "numeric"

# 解析 Q3 列名: "scan | roi | segment"
parts <- str_split_fixed(colnames(expr), " \\| ", 3)
meta <- data.frame(col=colnames(expr), scan=parts[,1],
                   roi=as.integer(parts[,2]), segment=parts[,3],
                   stringsAsFactors=FALSE)
seg2 <- seg %>% transmute(scan=`Scan name`, roi=as.integer(`ROI name`),
                          segment=`Segment name`, condition=as.character(Sample))
meta <- meta %>% left_join(seg2, by=c("scan","roi","segment"))

cat("== 未匹配到 condition 的 ROI 数:", sum(is.na(meta$condition)), "==\n")
stopifnot(all(!is.na(meta$condition)))   # join 必须全中
cat("== 分组样本量 ==\n"); print(table(meta$condition))

# 按时间排序 factor(control 为参照)
hrs <- suppressWarnings(as.numeric(str_extract(meta$condition, "\\d+")))
hrs[meta$condition=="control"] <- -1
lev <- unique(meta$condition[order(hrs)])
meta$condition <- factor(meta$condition, levels=lev)
meta$condition <- relevel(meta$condition, ref="control")

# ============ 2. limma 差异表达(各 CIP 时间点 vs control)============
logexpr <- log2(expr + 1)
design <- model.matrix(~0 + condition, data=meta)
colnames(design) <- levels(meta$condition)
fit <- lmFit(logexpr, design)

cips <- setdiff(levels(meta$condition), "control")
contr <- makeContrasts(contrasts=paste0(cips,"-control"), levels=design)
fit2 <- eBayes(contrasts.fit(fit, contr))

deg_list <- list(); deg_wb <- createWorkbook()
for(cc in cips){
  tt <- topTable(fit2, coef=paste0(cc,"-control"), number=Inf, sort.by="P")
  tt <- tt %>% rownames_to_column("gene")
  sig <- tt %>% filter(adj.P.Val < PADJ, abs(logFC) > LFC)
  deg_list[[cc]] <- sig$gene
  addWorksheet(deg_wb, cc); writeData(deg_wb, cc, tt)
  cat(sprintf("  %s: DEG=%d (up=%d, down=%d)\n", cc, nrow(sig),
              sum(sig$logFC>0), sum(sig$logFC<0)))
}
saveWorkbook(deg_wb, "results/DEG_limma_perTimepoint.xlsx", overwrite=TRUE)

# ============ 3. panel a: 时间点 DEG 韦恩 → 核心基因 ============
venn_sets <- if(length(cips)<=4) cips else intersect(c("CIP_3hr","CIP_12hr","CIP_24hr"), cips)
vlist <- deg_list[venn_sets]
cols <- c("#4E79A7","#F28E2B","#59A14F","#E15759")[seq_along(vlist)]
venn.diagram(vlist, filename="results/fig1_venn_timepoints.png",
             imagetype="png", height=1800, width=1800, resolution=300,
             fill=cols, alpha=0.55, cex=1.2, cat.cex=1.1, margin=0.12,
             main="AP time-course DEGs (vs control)")
core <- Reduce(intersect, vlist)
cat("== 核心交集 DEG 数:", length(core), "==\n")
if(length(core) < 10){ core <- unique(unlist(vlist)); cat("交集过小,改用并集, n=", length(core), "\n") }
writeLines(core, "results/core_DEG_genes.txt")

# ============ 4. panel b: GO/KEGG 富集(response to bacterium 在此)============
eg <- bitr(core, "SYMBOL", "ENTREZID", org.Mm.eg.db)
ego <- enrichGO(eg$ENTREZID, org.Mm.eg.db, ont="BP", keyType="ENTREZID",
                pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
write.csv(as.data.frame(ego), "results/GO_BP_enrich.csv", row.names=FALSE)

god <- as.data.frame(ego) %>% arrange(p.adjust) %>% head(20) %>%
  mutate(Description=factor(Description, levels=rev(Description)),
         neglogp=-log10(p.adjust),
         bacterium=grepl("bacter", Description, ignore.case=TRUE))
p_b <- ggplot(god, aes(neglogp, Description, fill=neglogp)) +
  geom_col() +
  geom_text(aes(label=ifelse(bacterium,"★","")), hjust=-0.3, color="red", size=5) +
  scale_fill_gradient(low="#FFE08A", high="#B2182B", name="-log10(padj)") +
  labs(x="-log10(p.adjust)", y=NULL, title="GO:BP enrichment (★ = 菌相关)") +
  theme_bw(base_size=12)
ggsave("results/fig2_GO_BP_bar.png", p_b, width=9, height=6, dpi=300)

bact <- as.data.frame(ego) %>% filter(grepl("bacter", Description, ignore.case=TRUE))
cat("== 命中的菌相关 GO 词条 ==\n"); print(bact[,c("ID","Description","p.adjust","Count")])

# ============ 5. panel c: KEGG 炎症通路 dotplot + 细胞因子 Sankey ============
kegg_ok <- tryCatch({
  ekegg <- enrichKEGG(eg$ENTREZID, organism="mmu", pvalueCutoff=0.1)
  ekegg <- setReadable(ekegg, org.Mm.eg.db, "ENTREZID"); TRUE
}, error=function(e){ cat("KEGG 联网失败(", conditionMessage(e), "),跳过 panel c;GO 部分已完成。\n"); FALSE })

if(kegg_ok){
  kdf <- as.data.frame(ekegg)
  write.csv(kdf, "results/KEGG_enrich.csv", row.names=FALSE)
  key <- "cytokine|chemokine|IL-17|TNF|Toll-like|NF-kappa|inflammat|NOD-like"
  kk <- kdf %>% filter(grepl(key, Description, ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(12) %>%
    mutate(M=as.numeric(sub("/.*","",BgRatio)),
           RichFactor=Count/M, neglogp=-log10(p.adjust),
           Description=factor(Description, levels=rev(Description)))
  p_c1 <- ggplot(kk, aes(RichFactor, Description, size=Count, color=neglogp)) +
    geom_point() +
    scale_color_gradient(low="#F6C445", high="#7A0177", name="-log10(padj)") +
    scale_size(range=c(3,10), name="Gene Number") +
    labs(x="Rich Factor", y=NULL, title="Inflammatory KEGG pathways") +
    theme_bw(base_size=12)
  ggsave("results/fig3_KEGG_dotplot.png", p_c1, width=9, height=6, dpi=300)

  # 细胞因子/趋化因子 gene→pathway Sankey
  fam <- "^(Ccl|Cxcl|Ccr|Cxcr|Il[0-9]|Tnf|Cx3c)"
  edges <- kk %>% transmute(pathway=as.character(Description), geneID) %>%
    tidyr::separate_rows(geneID, sep="/") %>%
    filter(grepl(fam, geneID)) %>% distinct()
  if(nrow(edges) > 3){
    p_c2 <- ggplot(edges, aes(axis1=geneID, axis2=pathway)) +
      geom_alluvium(aes(fill=pathway), width=1/12, alpha=0.7) +
      geom_stratum(width=1/12) +
      geom_text(stat="stratum", aes(label=after_stat(stratum)), size=2.8) +
      scale_x_discrete(limits=c("Gene","Pathway"), expand=c(.05,.05)) +
      theme_void(base_size=11) + theme(legend.position="none") +
      ggtitle("Cytokine/chemokine genes → KEGG pathways")
    ggsave("results/fig4_cytokine_sankey.png", p_c2, width=8, height=9, dpi=300)
  } else cat("细胞因子家族基因过少,跳过 Sankey。\n")
}
cat("\n==== 全部完成,结果在 results/ ====\n")
