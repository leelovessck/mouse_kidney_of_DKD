rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()

#打开必要的package
{
  if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(R.utils))install.packages("R.utils")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(stringr))install.packages("stringr")
  if(!require(enrichplot))install.packages("enrichplot")
  if(!require(msigdbr))install.packages("msigdbr")
  if(!require(GSVA))BiocManager::install("GSVA")
  if(!require(pheatmap))install.packages("pheatmap")
  if(!require(limma))BiocManager::install("limma")
  if(!require(BiocParallel))install.packages("BiocParallel")
  if(!require(ReactomePA))BiocManager::install("ReactomePA")
}

####修改下载协议#####################
#修改下载协议
R.utils::setOption("clusterProfiler.download.method","auto")



####载入数据###############
#载入数据
ALL.DEG <- read.csv("./result/5.3.2 WGCNA分析-EC/EC亚群hub基因-top25.csv")
rownames(ALL.DEG) <- ALL.DEG[,2]
colnames(ALL.DEG)[2] <- "SYMBOL"


DKD_rich <- ALL.DEG %>%     
  filter(module %in% c("M4", "M7", "M11", "M12", "M15", "M17", 
                       "M18", "M24", "M25", "M27")) %>%  
  group_by(module) %>%  
  arrange(desc(kME)) %>%  
  group_by(module, .add = TRUE) %>%  
  slice(1:100) %>%  
  ungroup()

CON_rich <- ALL.DEG %>%     
  filter(module %in% c("M1", "M3", "M10", "M14", "M16", 
                       "M19", "M20", "M22", "M23")) %>%  
  group_by(module) %>%  
  arrange(desc(kME)) %>%  
  group_by(module, .add = TRUE) %>%  
  slice(1:100) %>%  
  ungroup()


DKD_df <- bitr(DKD_rich$SYMBOL,
               fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
               OrgDb = org.Mm.eg.db)
DKD_df <- DKD_df %>% distinct(SYMBOL, .keep_all = T)

CON_df <- bitr(CON_rich$SYMBOL,
               fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
               OrgDb = org.Mm.eg.db)
CON_df <- CON_df %>% distinct(SYMBOL, .keep_all = T)



####KEGG富集分析（DKD）#################
#KEGG富集分析（DKD）
kegg <- enrichKEGG(unique(DKD_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.1 EC-DKD的KEGG.csv")



####GO（molecular function）（DKD）#################
#GO（molecular function）（DKD）
goMF <- enrichGO(DKD_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "1.2 EC-DKD的GO_MF.csv")



####GO（cell component）（DKD）######################
#GO（cell component）（DKD）
goCC <- enrichGO(DKD_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "1.3 EC-DKD的GO_CC.csv")



####GO（biological process）（DKD）##########################
#GO（biological process）（DKD）
goBP <- enrichGO(DKD_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "1.4 EC-DKD的GO_BP.csv")



####Reactome（DKD）##########################
#Reactome（DKD）
Reactome <- enrichPathway(unique(DKD_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "1.5 EC-DKD的Reactome.csv")



####KEGG富集分析（CON）#################
#KEGG富集分析（CON）
rm(kegg_category, kegg, goBP, goCC, goMF, Reactome)
kegg <- enrichKEGG(unique(CON_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.1 EC-CON的KEGG.csv")



####GO（molecular function）（CON）#################
#GO（molecular function）（CON）
goMF <- enrichGO(CON_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "2.2 EC-CON的GO_MF.csv")



####GO（cell component）（CON）######################
#GO（cell component）（CON）
goCC <- enrichGO(CON_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "2.3 EC-CON的GO_CC.csv")



####GO（biological process）（CON）##########################
#GO（biological process）（CON）
goBP <- enrichGO(CON_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "2.4 EC-CON的GO_BP.csv")



####Reactome（CON）##########################
#Reactome（CON）
Reactome <- enrichPathway(unique(CON_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.5 EC-CON的Reactome.csv")
