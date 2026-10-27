suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(stringr)
  library(limma); library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2)
  library(ggalluvial); library(tidyr)
})
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")

# 复用已存的 12hr DEG 结果(上一脚本已写入 xlsx)
tt <- read_excel("results/DEG_limma_perTimepoint.xlsx", sheet="CIP_12hr")
sig <- tt %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("12hr DEG 数:", nrow(sig), "\n")
eg <- bitr(sig$gene, "SYMBOL", "ENTREZID", org.Mm.eg.db)

ekegg <- enrichKEGG(eg$ENTREZID, organism="mmu", pvalueCutoff=0.1)
ekegg <- setReadable(ekegg, org.Mm.eg.db, "ENTREZID")
kdf <- as.data.frame(ekegg)
write.csv(kdf, "results/KEGG_12hr.csv", row.names=FALSE)
cat("KEGG 富集通路总数:", nrow(kdf), "\n")
cat("== 前 20 条通路名(确认有没有炎症相关)==\n"); print(head(kdf$Description, 20))

# 宽松匹配炎症/免疫通路;匹配不到就退回 top15 全部
key <- "cytokine|chemokine|IL-17|IL17|TNF|Toll|NF-kappa|NFkB|inflammat|NOD|Th17|complement|phagosome|Leukocyte|Chemokine"
kk <- kdf %>% filter(grepl(key, Description, ignore.case=TRUE))
if(nrow(kk) < 3){
  cat("!! 炎症关键词命中<3,退回展示 top15 全部通路\n")
  kk <- kdf %>% arrange(p.adjust) %>% head(15)
} else kk <- kk %>% arrange(p.adjust) %>% head(15)

kk <- kk %>% mutate(M=as.numeric(sub("/.*","",BgRatio)),
                    RichFactor=Count/M, neglogp=-log10(p.adjust),
                    Description=factor(Description, levels=rev(Description)))
p_c1 <- ggplot(kk, aes(RichFactor, Description, size=Count, color=neglogp)) +
  geom_point() +
  scale_color_gradient(low="#F6C445", high="#7A0177", name="-log10(padj)") +
  scale_size(range=c(3,10), name="Gene Number") +
  labs(x="Rich Factor", y=NULL, title="KEGG pathways @ CIP_12hr (inflammation peak)") +
  theme_bw(base_size=12) + theme(axis.text.y=element_text(size=10))
ggsave("results/fig3_KEGG_dotplot.png", p_c1, width=10, height=6.5, dpi=300)

# Sankey:趋化/细胞因子基因 → 通路
fam <- "^(Ccl|Cxcl|Ccr|Cxcr|Il[0-9]|Tnf|Cx3c|Nfkb|Tlr)"
edges <- kk %>% transmute(pathway=as.character(Description), geneID) %>%
  separate_rows(geneID, sep="/") %>% filter(grepl(fam, geneID)) %>% distinct()
cat("Sankey 边数:", nrow(edges), "\n")
if(nrow(edges) >= 3){
  p_c2 <- ggplot(edges, aes(axis1=geneID, axis2=pathway)) +
    geom_alluvium(aes(fill=pathway), width=1/12, alpha=0.7) +
    geom_stratum(width=1/12) +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=2.8) +
    scale_x_discrete(limits=c("Gene","Pathway"), expand=c(.05,.05)) +
    theme_void(base_size=11) + theme(legend.position="none") +
    ggtitle("Cytokine/chemokine genes -> KEGG pathways (12hr)")
  ggsave("results/fig4_cytokine_sankey.png", p_c2, width=8, height=9, dpi=300)
} else cat("边仍过少,Sankey 跳过(用 dotplot 即可)\n")
cat("==== 完成 ====\n")
