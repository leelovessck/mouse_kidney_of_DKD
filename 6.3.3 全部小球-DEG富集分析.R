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
ALL.DEG <- read.csv("./result/5.2 全部细胞DEG/全部细胞的DEG.csv")
rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"



####全部DEG####################
#全部DEG
ALL.DEG_rich <- ALL.DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

ALL.DEG_df <- bitr(ALL.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Mm.eg.db)
ALL.DEG_df <- ALL.DEG_df %>% distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-ALL#################
#KEGG富集分析-ALL
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
write.csv(kegg@result, "1.1 全部细胞的KEGG（all）.csv")



####GO（molecular function）-ALL#################
#GO（molecular function）-ALL
goMF <- enrichGO(ALL.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "1.2 全部细胞的GO_MF（all）.csv")



####GO（cell component）-ALL######################
#GO（cell component）-ALL
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
write.csv(goCC@result, "1.3 全部细胞的GO_CC（all）.csv")



####GO（biological process）-ALL##########################
#GO（biological process）-ALL
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
write.csv(goBP@result, "1.4 全部细胞的GO_BP（all）.csv")



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
write.csv(Reactome@result, "1.5 全部细胞的Reactome（all）.csv")



####上调DEG####################
#上调DEG
rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)

UP.DEG_rich <- ALL.DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

UP.DEG_df <- bitr(UP.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
UP.DEG_df <- UP.DEG_df %>% distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-上调#################
#KEGG富集分析-上调
kegg <- enrichKEGG(unique(UP.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.1 全部细胞的KEGG（up）.csv")



####GO（molecular function）-上调#################
#GO（molecular function）-上调
goMF <- enrichGO(UP.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "2.2全部细胞的GO_MF（up）.csv")



####GO（cell component）-上调######################
#GO（cell component）-上调
goCC <- enrichGO(UP.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "2.3 全部细胞的GO_CC（up）.csv")



####GO（biological process）-上调##########################
#GO（biological process）-上调
goBP <- enrichGO(UP.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "2.4 全部细胞的GO_BP（up）.csv")



####Reactome-上调##########################
#Reactome-上调
Reactome <- enrichPathway(unique(UP.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.5 全部细胞的Reactome（up）.csv")



####下调DEG####################
#下调DEG
rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)

DOWN.DEG_rich <- ALL.DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

DOWN.DEG_df <- bitr(DOWN.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                    OrgDb = org.Mm.eg.db)
DOWN.DEG_df <- DOWN.DEG_df %>% distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-下调#################
#KEGG富集分析-下调
kegg <- enrichKEGG(unique(DOWN.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "3.1 全部细胞的KEGG（down）.csv")



####GO（molecular function）-下调#################
#GO（molecular function）-下调
goMF <- enrichGO(DOWN.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "3.2 全部细胞的GO_MF（down）.csv")



####GO（cell component）-下调######################
#GO（cell component）-下调
goCC <- enrichGO(DOWN.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "3.3 全部细胞的GO_CC（down）.csv")



####GO（biological process）-下调##########################
#GO（biological process）-下调
goBP <- enrichGO(DOWN.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "3.4 全部细胞的GO_BP（down）.csv")



####Reactome-下调##########################
#Reactome-下调
Reactome <- enrichPathway(unique(DOWN.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "3.5 全部细胞的Reactome（down）.csv")
