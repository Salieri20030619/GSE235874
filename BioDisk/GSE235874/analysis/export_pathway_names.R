suppressPackageStartupMessages({ library(dplyr); library(stringr) })
setwd("/mnt/BioDisk/wangmingyuan/GSE235874")

# 复现图里的 abbr() 缩写函数(与主脚本完全一致)
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

# ---- C: GO:BP 核心 DEG,取图上展示的 top15 ----
go <- read.csv("results/GO_BP_core.csv", stringsAsFactors=FALSE) %>%
  arrange(p.adjust) %>% head(15) %>%
  transmute(panel="C (GO:BP)", 缩写=abbr(Description), 全称=Description,
            ID=ID, p.adjust=signif(p.adjust,3), Count=Count)

# ---- D: KEGG @12hr,复现图上的过滤+top14 ----
kegg <- read.csv("results/KEGG_12hr.csv", stringsAsFactors=FALSE)
noise <- "Ribosome|Coronavirus|Estrogen|Alcoholism|Oxytocin|Proteoglycans in cancer|Protein digestion|Pancreatic secretion"
ikey  <- "bact|infection|Salmonella|trap|Toll|NOD|phagos|Leukocyte|Antigen|Neutrophil|Complement|Platelet|adhesion|cytoskeleton|junction|migration|Autophagy|Chemokine|cytokine|inflammat"
kk <- kegg %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>%
  filter(grepl(ikey,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
if(nrow(kk)<6) kk <- kegg %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
kegg_out <- kk %>% transmute(panel="D (KEGG)", 缩写=abbr(Description), 全称=Description,
            ID=ID, p.adjust=signif(p.adjust,3), Count=Count)

both <- bind_rows(go, kegg_out)
write.csv(both, "results/pathway_abbr_fullname.csv", row.names=FALSE)
cat("==== 导出", nrow(both), "条通路对照 ====\n")
cat("C (GO):", nrow(go), "条  |  D (KEGG):", nrow(kegg_out), "条\n\n")
print(both[, c("panel","缩写","全称")], right=FALSE)
