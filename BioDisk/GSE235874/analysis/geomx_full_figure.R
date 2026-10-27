# ============================================================
# GSE235874 GeoMx WTA CIP 急性胰腺炎 9-panel(统一自定义配色)
# ============================================================
suppressPackageStartupMessages({
  need_cran <- c("readxl","dplyr","tibble","stringr","tidyr","ggplot2","ggVennDiagram",
                 "patchwork","pheatmap","ggplotify","ggrepel","RColorBrewer","cowplot")
  need_bioc <- c("limma","clusterProfiler","org.Mm.eg.db","enrichplot","GSVA")
  mc <- need_cran[!need_cran %in% rownames(installed.packages())]
  mb <- need_bioc[!need_bioc %in% rownames(installed.packages())]
  if(length(mc)) install.packages(mc, repos="https://cloud.r-project.org")
  if(length(mb)){ if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
                  BiocManager::install(mb, update=FALSE, ask=FALSE) }
  library(readxl); library(dplyr); library(tibble); library(stringr); library(tidyr)
  library(ggplot2); library(ggVennDiagram); library(patchwork); library(pheatmap)
  library(ggplotify); library(ggrepel); library(RColorBrewer); library(cowplot)
  library(limma); library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot); library(GSVA)
})
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")
dir.create("results/panels", showWarnings=FALSE, recursive=TRUE)
set.seed(1)
BASE <- 13; TITLE <- 14; STAR <- "*"

# ========== 统一自定义色板 ==========
PAL_BLUE <- "#00C5FD"; PAL_PURPLE <- "#E16CFC"; PAL_GREEN <- "#51E6AD"
PAL_MAGENTA <- "#EE3BA5"; PAL_RED <- "#D83010"; PAL_ORANGE <- "#FFC565"
PAL_PINK <- "#FF9EC4"; PAL_NAVY <- "#1F4E9C"
PAL2 <- c(PAL_BLUE, PAL_PURPLE)                                              # 2 组
PAL3 <- c(PAL_BLUE, PAL_PURPLE, PAL_GREEN)                                   # 3 组
PAL6 <- c(PAL_MAGENTA, PAL_RED, PAL_ORANGE, PAL_PURPLE, PAL_BLUE, PAL_GREEN) # 6 组
PAL8 <- c(PAL6, PAL_PINK, PAL_NAVY)                                          # 8 组

bold_theme <- theme(
  text=element_text(face="bold"), axis.text=element_text(face="bold", color="black"),
  axis.title=element_text(face="bold"), legend.text=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.title=element_text(face="bold", size=TITLE))

abbr <- function(x){
  x <- as.character(x)
  reps <- c("antimicrobial"="antimicr.","immune response"="imm.resp.","response"="resp.",
            "regulation"="reg.","negative"="neg.","positive"="pos.","organization"="org.",
            "supramolecular fiber"="supramol.fiber","protein-containing complex"="prot.complex",
            "defense"="def.","bacterium"="bact.","Gram-positive"="Gram+","Gram-negative"="Gram-",
            "processing and presentation"="proc.&present.","transendothelial migration"="transendo.migr.",
            "differentiation"="diff.","development"="dev.","signaling pathway"="signaling",
            "involved in"="in","targeting to vacuole"="to vacuole","humoral"="hum.")
  for(k in names(reps)) x <- gsub(k, reps[[k]], x, ignore.case=TRUE)
  str_trunc(x, 32)
}

LFC <- 1.0; PADJ <- 0.05
q3_file  <- "GSE235874_QC_BQC_Q3Norm_matrixonly.xlsx"
seg_file <- "GSE235874_segmentSummary_readTypes.xlsx"
retry <- function(expr,times=5,wait=6,what="op"){for(i in seq_len(times)){
  r<-tryCatch(force(expr),error=function(e){message(sprintf("[%s]%d失败:%s",what,i,conditionMessage(e)));NULL})
  if(!is.null(r))return(r); if(i<times){message(sprintf("[%s]%ds后重试",what,wait));Sys.sleep(wait)}}; NULL}

# ============ 0. 读入 ============
q3 <- read_excel(q3_file); seg <- read_excel(seg_file)
expr <- as.data.frame(q3); rownames(expr) <- make.unique(as.character(expr$TargetName))
expr$TargetName <- NULL; expr <- as.matrix(expr); storage.mode(expr) <- "numeric"
logexpr <- log2(expr + 1)
parts <- str_split_fixed(colnames(expr), " \\| ", 3)
meta <- data.frame(col=colnames(expr), scan=parts[,1], roi=as.integer(parts[,2]),
                   segment=parts[,3], stringsAsFactors=FALSE)
seg2 <- seg %>% transmute(scan=`Scan name`, roi=as.integer(`ROI name`),
                          segment=`Segment name`, condition=as.character(Sample))
meta <- meta %>% left_join(seg2, by=c("scan","roi","segment")); stopifnot(all(!is.na(meta$condition)))
parse_hr <- function(x){ if(x=="control") return(0); n<-as.numeric(str_extract(x,"\\d+"))
                         if(grepl("hr",x)) n else if(grepl("d",x)) n*24 else NA }
tp_levels <- c("control","CIP_3hr","CIP_12hr","CIP_24hr","CIP_48hr","CIP_72hr","CIP_5d","CIP_7d","CIP_14d","CIP_30d")
tp_levels <- tp_levels[tp_levels %in% unique(meta$condition)]
meta$condition <- factor(meta$condition, levels=tp_levels); meta <- meta[order(meta$condition),]
logexpr <- logexpr[, meta$col]; expr <- expr[, meta$col]; hr_map <- sapply(tp_levels, parse_hr)

# ============ 1. limma ============
design <- model.matrix(~0 + condition, data=meta); colnames(design) <- levels(meta$condition)
fit <- lmFit(logexpr, design); cips <- setdiff(levels(meta$condition), "control")
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=paste0(cips,"-control"), levels=design)))
deg_tab <- list(); deg_genes <- list()
for(cc in cips){
  tt <- topTable(fit2, coef=paste0(cc,"-control"), number=Inf, sort.by="P") %>% rownames_to_column("gene")
  deg_tab[[cc]] <- tt; deg_genes[[cc]] <- tt %>% filter(adj.P.Val<PADJ, abs(logFC)>LFC) %>% pull(gene)
}

# ============ A (2组: 蓝紫) ============
degN <- lapply(cips, function(cc){ tt <- deg_tab[[cc]] %>% filter(adj.P.Val<PADJ, abs(logFC)>LFC)
  data.frame(tp=cc, up=sum(tt$logFC>0), down=sum(tt$logFC<0))}) %>% bind_rows() %>%
  mutate(tp=factor(tp,levels=cips), down=-down) %>% pivot_longer(c(up,down),names_to="dir",values_to="n")
pA <- ggplot(degN, aes(tp,n,fill=dir)) + geom_col(width=0.7) + geom_hline(yintercept=0,color="grey40") +
  scale_fill_manual(values=c(up=PAL_PURPLE, down=PAL_BLUE),name="") +
  labs(x=NULL,y="DEG count",title="A. DEG dynamics across AP time course") +
  theme_bw(base_size=BASE) + bold_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,face="bold"),
        legend.position="top", legend.direction="horizontal")

# ============ B (韦恩渐变, 不动) ============
vsets <- intersect(c("CIP_3hr","CIP_12hr","CIP_24hr"), cips); vlist <- deg_genes[vsets]; names(vlist) <- vsets
pB <- ggVennDiagram(vlist, label="count", label_alpha=0, set_size=5, label_size=5, edge_size=0.6) +
  scale_fill_gradient(low="#F0F4F8", high=PAL_ORANGE, guide="none") +
  scale_color_manual(values=c(PAL_BLUE, PAL_ORANGE, PAL_GREEN)) +
  coord_cartesian(clip="off") + labs(title="B. Core DEGs (acute phase)") +
  theme(plot.title=element_text(face="bold", size=TITLE, hjust=0.5), plot.margin=margin(15,25,15,25))
core <- Reduce(intersect, vlist); writeLines(core, "results/core_DEG_genes.txt")

# ============ C (2组: 蓝紫) ============
eg <- bitr(core,"SYMBOL","ENTREZID",org.Mm.eg.db)
ego <- enrichGO(eg$ENTREZID, org.Mm.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
god <- as.data.frame(ego) %>% arrange(p.adjust) %>% head(15) %>%
  mutate(Description=str_wrap(Description,width=30), Description=factor(Description,levels=rev(Description)),
         neglogp=-log10(p.adjust), bact=grepl("bact|antimicr",Description,ignore.case=TRUE))
pC <- ggplot(god, aes(neglogp,Description,fill=bact)) + geom_col() +
  scale_fill_manual(values=c(`TRUE`=PAL_PURPLE,`FALSE`=PAL_BLUE),labels=c("other","microbe-rel."),name="") +
  labs(x="-log10(p.adjust)",y=NULL,title="C. GO:BP of core DEGs") +
  theme_bw(base_size=BASE) + bold_theme +
  theme(axis.text.y=element_text(size=BASE-2,face="bold"), legend.position="top", legend.direction="horizontal")

# ============ D (连续色阶, 不动) ============
eg12 <- bitr(deg_genes[["CIP_12hr"]],"SYMBOL","ENTREZID",org.Mm.eg.db)
invisible(retry(clusterProfiler::download_KEGG("mmu"), times=5, wait=6, what="KEGG预取"))
ekegg <- retry(setReadable(enrichKEGG(eg12$ENTREZID,"mmu",pvalueCutoff=0.15), org.Mm.eg.db,"ENTREZID"),
               times=5, wait=8, what="enrichKEGG")
if(!is.null(ekegg) && nrow(as.data.frame(ekegg))>0){
  kdf <- as.data.frame(ekegg); write.csv(kdf,"results/KEGG_12hr.csv",row.names=FALSE)
  noise <- "Ribosome|Coronavirus|Estrogen|Alcoholism|Oxytocin|Proteoglycans in cancer|Protein digestion|Pancreatic secretion"
  ikey <- "bact|infection|Salmonella|trap|Toll|NOD|phagos|Leukocyte|Antigen|Neutrophil|Complement|Platelet|adhesion|cytoskeleton|junction|migration|Autophagy|Chemokine|cytokine|inflammat"
  pathkey <- "bact|Salmonella|infection|trap|phagos|virus|herpes|immunodeficiency|cytomegalo|NOD|Toll"
  kk <- kdf %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>%
    filter(grepl(ikey,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
  if(nrow(kk)<6) kk <- kdf %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
  kk <- kk %>% mutate(M=as.numeric(sub("/.*","",BgRatio)),RichFactor=Count/M,neglogp=-log10(p.adjust),
           pathogen=grepl(pathkey,Description,ignore.case=TRUE),
           Description=str_wrap(Description,width=30),Description=factor(Description,levels=rev(Description)))
  pD <- ggplot(kk, aes(RichFactor,Description,size=Count,color=neglogp)) + geom_point() +
    geom_text(aes(label=ifelse(pathogen,STAR,"")),color="red",size=6,fontface="bold",nudge_x=0.003) +
    scale_color_gradient(low=PAL_ORANGE,high=PAL_PURPLE,name="-log10\n(padj)") +
    scale_size(range=c(2.5,8),name="Count") +
    labs(x="Rich Factor",y=NULL,title="D. KEGG @12hr (* = pathogen-related)") +
    theme_bw(base_size=BASE) + bold_theme + theme(axis.text.y=element_text(size=BASE-4,face="bold"),
      legend.key.height=unit(11,"pt"))
} else pD <- ggplot()+annotate("text",1,1,label="KEGG unavailable",size=5,fontface="bold")+theme_void()+ggtitle("D. KEGG @12hr")

# ============ GO gene sets ============
go_terms <- c("defense response to bacterium"="GO:0042742","antimicrobial humoral response"="GO:0019730",
              "acute inflammatory response"="GO:0002526","neutrophil chemotaxis"="GO:0030593",
              "response to wounding"="GO:0009611","leukocyte migration"="GO:0050900")
gene_sets <- lapply(go_terms, function(gid){
  eids <- tryCatch(AnnotationDbi::select(org.Mm.eg.db,keys=gid,keytype="GOALL",columns="SYMBOL")$SYMBOL,
                   error=function(e) character(0)); intersect(unique(eids), rownames(logexpr))})
gene_sets <- gene_sets[sapply(gene_sets,length)>=3]

# ============ E (6组: 六色板, 图例2列大字) ============
gp <- try(GSVA::gsvaParam(logexpr, gene_sets, kcdf="Gaussian"), silent=TRUE)
gsva_mat <- if(!inherits(gp,"try-error")) GSVA::gsva(gp) else GSVA::gsva(logexpr,gene_sets,kcdf="Gaussian",verbose=FALSE)
gsva_long <- as.data.frame(t(gsva_mat)) %>% rownames_to_column("col") %>%
  left_join(meta[,c("col","condition")],by="col") %>%
  pivot_longer(-c(col,condition),names_to="pathway",values_to="score") %>%
  group_by(condition,pathway) %>% summarise(mean=mean(score),se=sd(score)/sqrt(n()),.groups="drop") %>%
  mutate(tp_idx=as.integer(factor(condition,levels=tp_levels)))
yr <- range(gsva_long$mean-gsva_long$se, gsva_long$mean+gsva_long$se, na.rm=TRUE)
pE <- ggplot(gsva_long, aes(tp_idx,mean,color=pathway)) +
  geom_ribbon(aes(ymin=mean-se,ymax=mean+se,fill=pathway),alpha=0.10,color=NA) +
  geom_line(linewidth=1.2)+geom_point(size=2.4) +
  scale_x_continuous(breaks=seq_along(tp_levels),labels=tp_levels, expand=expansion(mult=0.02)) +
  scale_y_continuous(limits=yr, expand=expansion(mult=0.03)) +
  scale_color_manual(values=PAL6)+scale_fill_manual(values=PAL6) +
  labs(x=NULL,y="GSVA score",title="E. Immune/antimicrobial activity over time") +
  guides(color=guide_legend(ncol=2,byrow=TRUE),fill=guide_legend(ncol=2,byrow=TRUE)) +
  theme_bw(base_size=BASE) + bold_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=BASE-3,face="bold"),
    legend.position="bottom", legend.text=element_text(size=BASE-1,face="bold"),
    legend.key.height=unit(8,"pt"), legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"),
    legend.title=element_text(size=BASE-1,face="bold"))

# ============ F (热图 RdBu, 不动) ============
foc_genes <- unique(unlist(gene_sets[c("defense response to bacterium","antimicrobial humoral response","acute inflammatory response")]))
foc_genes <- intersect(foc_genes, unique(unlist(deg_genes)))
tp_mean <- sapply(tp_levels, function(tp) rowMeans(logexpr[,meta$col[meta$condition==tp],drop=FALSE]))
hm <- tp_mean[foc_genes,,drop=FALSE]; hm <- hm[order(-apply(hm,1,var)),,drop=FALSE]; hm <- head(hm,40)
ph <- pheatmap(t(scale(t(hm))), cluster_cols=FALSE, color=colorRampPalette(c("#00C5FD","#FFFFFF","#E16CFC"))(50),
         fontsize_row=9, fontsize_col=11, border_color=NA, angle_col=45,
         main="F. Antimicrobial/inflammatory genes over time", silent=TRUE)
pF <- as.ggplot(ph)

# ============ G (2组+ns: 蓝紫灰) ============
vol <- deg_tab[["CIP_12hr"]] %>% mutate(sig=ifelse(adj.P.Val<PADJ & abs(logFC)>LFC, ifelse(logFC>0,"up","down"),"ns"))
lab <- vol %>% filter(gene %in% foc_genes, sig!="ns") %>% arrange(adj.P.Val) %>% head(15)
pG <- ggplot(vol, aes(logFC,-log10(adj.P.Val),color=sig)) + geom_point(size=2,alpha=0.6) +
  geom_vline(xintercept=c(-LFC,LFC),linetype=2,color="grey60") + geom_hline(yintercept=-log10(PADJ),linetype=2,color="grey60") +
  geom_text_repel(data=lab,aes(label=gene),size=3.4,fontface="bold",max.overlaps=20,color="black") +
  scale_color_manual(values=c(up=PAL_PURPLE,down=PAL_BLUE,ns="grey80")) +
  annotate("text", x=max(vol$logFC)*0.7, y=max(-log10(vol$adj.P.Val))*0.95, label="up", color=PAL_PURPLE, fontface="bold", size=5) +
  annotate("text", x=min(vol$logFC)*0.7, y=max(-log10(vol$adj.P.Val))*0.95, label="down", color=PAL_BLUE, fontface="bold", size=5) +
  labs(title="G. Volcano @12hr (peak)",x="log2FC",y="-log10(padj)") +
  theme_bw(base_size=BASE) + bold_theme + theme(legend.position="none")

# ============ H (Reds 渐变, 不动) ============
bac_go <- as.data.frame(ego) %>% filter(grepl("bact|antimicr",Description,ignore.case=TRUE)) %>%
  arrange(p.adjust) %>% head(5)
if(nrow(bac_go)>0){
  hedges <- bac_go %>% transmute(GO=str_wrap(Description,width=28), geneID) %>%
    separate_rows(geneID,sep="/") %>% distinct() %>% mutate(present=1)
  keep_g <- hedges %>% count(geneID) %>% arrange(-n) %>% head(25) %>% pull(geneID)
  hedges <- hedges %>% filter(geneID %in% keep_g)
  g_order <- hedges %>% count(geneID) %>% arrange(n) %>% pull(geneID)
  pH <- ggplot(hedges, aes(x=factor(geneID,levels=g_order), y=GO)) +
    geom_tile(aes(fill=GO),color="white",show.legend=FALSE) + scale_fill_manual(values=colorRampPalette(c("#F3E6FC","#E16CFC","#8E24AA"))(5)) +
    labs(x=NULL,y=NULL,title="H. Microbe-response GO membership") +
    theme_minimal(base_size=BASE) + bold_theme +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=BASE-4,face="bold"),
          axis.text.y=element_text(size=BASE-2,face="bold"), panel.grid=element_blank())
} else pH <- ggplot()+annotate("text",1,1,label="no microbe GO",fontface="bold")+theme_void()+ggtitle("H.")

# ============ I (8组: 八色板) ============
key_genes <- intersect(c("Reg3a","Reg3b","Reg3g","Reg2","Lyz2","Dmbt1","Lcn2","Hamp2"), rownames(logexpr))
dfI <- logexpr[key_genes,,drop=FALSE] %>% t() %>% as.data.frame() %>% rownames_to_column("col") %>%
  left_join(meta[,c("col","condition")], by="col") %>%
  pivot_longer(-c(col,condition), names_to="gene", values_to="expr") %>%
  group_by(condition, gene) %>% summarise(mean=mean(expr), se=sd(expr)/sqrt(n()), .groups="drop") %>%
  mutate(tp_idx=as.integer(factor(condition, levels=tp_levels)))
pI <- ggplot(dfI, aes(tp_idx, mean, color=gene)) +
  geom_line(linewidth=1.2) + geom_point(size=2.2) +
  scale_x_continuous(breaks=seq_along(tp_levels), labels=tp_levels, expand=expansion(mult=0.02)) +
  scale_color_manual(values=PAL8) +
  labs(x=NULL, y="log2 expression", title="I. Antimicrobial peptide genes over time") +
  guides(color=guide_legend(ncol=2)) +
  theme_bw(base_size=BASE) + bold_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=BASE-3,face="bold"),
        legend.position="bottom", legend.text=element_text(size=BASE-4,face="bold"),
        legend.key.height=unit(9,"pt"), legend.title=element_blank(),
        legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"))

# ============ 拼图 ============
dir.create("results/tiff_panels", showWarnings=FALSE)
.wp <- list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF,G=pG,H=pH,I=pI)
for(.nm in names(.wp)){
  tiff(sprintf("results/tiff_panels/WTA_panel_%s.tiff",.nm), width=7,height=6,units="in",res=600,compression="lzw")
  print(.wp[[.nm]]); dev.off()
}
cat(">> WTA 9 panel TIFF 完成\n")
big <- cowplot::plot_grid(pA,pB,pC,pD,pE,pF,pG,pH,pI, ncol=3, align="none")
title_gg <- ggplot()+labs(title="GSE235874 GeoMx WTA - Cerulein-induced acute pancreatitis: onset, progression & antimicrobial immunity")+theme_void()+theme(plot.title=element_text(face="bold",size=18,hjust=0.01))
big <- cowplot::plot_grid(title_gg, big, ncol=1, rel_heights=c(0.03,1))
big <- cowplot::ggdraw() + cowplot::draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA))) + cowplot::draw_plot(big)
ggsave("results/Figure_GSE235874_9panel_600dpi.png", big, width=22, height=22, dpi=600, limitsize=FALSE)
ggsave("results/Figure_GSE235874_9panel.pdf", big, width=22, height=22, limitsize=FALSE)
cat("\n==== 完成 ====\n核心DEG:",length(core)," 菌相关GO:",nrow(bac_go),
    " GSVA基因集:",length(gene_sets)," I基因数:",length(key_genes),
    " KEGG:",ifelse(is.null(ekegg),"失败","成功"),"\n")
