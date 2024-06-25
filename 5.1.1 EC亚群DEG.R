rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")



####全部EC############################
#全部EC
EC <- subset(mouse_all, usetype == "EC")
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
table(EC$usetype)
EC_DEG <- FindMarkers(EC, 
                      min.pct = 0.10, 
                      logfc.threshold = 0.10,
                      group.by = "group",
                      ident.1 = "DKD",
                      ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(EC_DEG, "全部EC的DEG.csv")

EnhancedVolcano(EC_DEG,
                lab = rownames(EC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "全部EC的DEG（2倍）")
  ##得到1.全部EC的DEG（2倍）

EnhancedVolcano(EC_DEG,
                lab = rownames(EC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                title = "全部EC的DEG（4倍）")
  ##得到2.全部EC的DEG（4倍）
