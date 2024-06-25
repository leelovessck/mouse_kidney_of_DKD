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
ALL.DEG <- read.csv("./result/5.3.1 WGCNA分析/hub基因-all.csv")
rownames(ALL.DEG) <- ALL.DEG[,2]
colnames(ALL.DEG)[2] <- "SYMBOL"



####全部DEG####################
#全部DEG
ALL.DEG_rich <- ALL.DEG %>%     
  filter(module %in% c("M1", "M8")) %>%  
  group_by(module) %>%  
  arrange(desc(kME)) %>%  
  group_by(module, .add = TRUE) %>%  
  slice(1:100) %>%  
  ungroup() 


ALL.DEG_df <- bitr(ALL.DEG_rich$SYMBOL,
                   fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Mm.eg.db)
ALL.DEG_df <- ALL.DEG_df %>% distinct(SYMBOL, .keep_all = T)



####KEGG富集分析#################
#KEGG富集分析
kegg <- enrichKEGG(unique(ALL.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.EC模块的KEGG.csv")



####GO（molecular function）#################
#GO（molecular function）
goMF <- enrichGO(ALL.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "2.EC模块的GO_MF.csv")



####GO（cell component）######################
#GO（cell component）
goCC <- enrichGO(ALL.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "3.EC模块的GO_CC.csv")



####GO（biological process）##########################
#GO（biological process）
goBP <- enrichGO(ALL.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "4.EC模块的GO_BP.csv")



####Reactome-ALL##########################
#Reactome-ALL
Reactome <- enrichPathway(unique(ALL.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "5.EC模块的Reactome（all）.csv")
