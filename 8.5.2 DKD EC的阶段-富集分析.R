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
all_marker <- read.csv("./result/8.5.1 DKD EC的阶段-基因/EC时序分组的marker.csv")
all_DEG <- read.csv("./result/8.5.1 DKD EC的阶段-基因/EC时序分组的DEG.csv")



####C1的marker####################
#C1的marker
C1.marker_rich <- all_marker %>%     
  filter(cluster == "1") %>%
  filter(p_val < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%  
  head(500) %>% 
  select(gene) %>%
  unlist() 

C1.DEG_df <- bitr(C1.marker_rich,
                  fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
C1.DEG_df <- C1.DEG_df %>% 
  distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-C1#################
#KEGG富集分析-C1
kegg <- enrichKEGG(unique(C1.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.1 时序基因的KEGG（C1）.csv")



####GO（molecular function）-C1#################
#GO（molecular function）-C1
goMF <- enrichGO(C1.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "1.2 时序基因的GO_MF（C1）.csv")



####GO（cell component）-C1######################
#GO（cell component）-C1
goCC <- enrichGO(C1.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "1.3 时序基因的GO_CC（C1）.csv")



####GO（biological process）-C1##########################
#GO（biological process）-C1
goBP <- enrichGO(C1.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "1.4 时序基因的GO_BP（C1）.csv")



####Reactome-C1##########################
#Reactome-C1
Reactome <- enrichPathway(unique(C1.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "1.5 时序基因的Reactome（C1）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)



####C2的marker####################
#C2的marker
C2.marker_rich <- all_marker %>%     
  filter(cluster == "2") %>%
  filter(p_val < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%  
  head(500) %>% 
  select(gene) %>%
  unlist() 

C2.DEG_df <- bitr(C2.marker_rich,
                  fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
C2.DEG_df <- C2.DEG_df %>% 
  distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-C2#################
#KEGG富集分析-C2
kegg <- enrichKEGG(unique(C2.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.1 时序基因的KEGG（C2）.csv")



####GO（molecular function）-C2#################
#GO（molecular function）-C2
goMF <- enrichGO(C2.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "2.2 时序基因的GO_MF（C2）.csv")



####GO（cell component）-C2######################
#GO（cell component）-C2
goCC <- enrichGO(C2.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "2.3 时序基因的GO_CC（C2）.csv")



####GO（biological process）-C2##########################
#GO（biological process）-C2
goBP <- enrichGO(C2.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "2.4 时序基因的GO_BP（C2）.csv")



####Reactome-C2##########################
#Reactome-C2
Reactome <- enrichPathway(unique(C2.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.5 时序基因的Reactome（C2）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)



####C3的marker####################
#C3的marker
C3.marker_rich <- all_marker %>%     
  filter(cluster == "2") %>%
  filter(p_val < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%  
  head(500) %>% 
  select(gene) %>%
  unlist() 

C3.DEG_df <- bitr(C3.marker_rich,
                  fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
C3.DEG_df <- C3.DEG_df %>% 
  distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-C3#################
#KEGG富集分析-C3
kegg <- enrichKEGG(unique(C3.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "3.1 时序基因的KEGG（C3）.csv")



####GO（molecular function）-C3#################
#GO（molecular function）-C3
goMF <- enrichGO(C3.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "3.2 时序基因的GO_MF（C3）.csv")



####GO（cell component）-C3######################
#GO（cell component）-C3
goCC <- enrichGO(C3.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "3.3 时序基因的GO_CC（C3）.csv")



####GO（biological process）-C3##########################
#GO（biological process）-C3
goBP <- enrichGO(C3.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "3.4 时序基因的GO_BP（C3）.csv")



####Reactome-C3##########################
#Reactome-C3
Reactome <- enrichPathway(unique(C3.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "3.5 时序基因的Reactome（C3）.csv")



####全部DEG####################
#全部DEG
ALL.DEG_rich <- all_DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  select(X) %>%
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
write.csv(kegg@result, "4.1 时序EC的KEGG（all）.csv")



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
write.csv(goMF@result, "4.2 时序EC的GO_MF（all）.csv")



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
write.csv(goCC@result, "4.3 时序EC的GO_CC（all）.csv")



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
write.csv(goBP@result, "4.4 时序EC的GO_BP（all）.csv")



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
write.csv(Reactome@result, "4.5 时序EC的Reactome（all）.csv")



####上调DEG####################
#上调DEG
rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)

UP.DEG_rich <- all_DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  select(X) %>%
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
write.csv(kegg@result, "5.1 时序EC的KEGG（up）.csv")



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
write.csv(goMF@result, "5.2 时序EC的GO_MF（up）.csv")



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
write.csv(goCC@result, "5.3 时序EC的GO_CC（up）.csv")



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
write.csv(goBP@result, "5.4 时序EC的GO_BP（up）.csv")



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
write.csv(Reactome@result, "5.5 时序EC的Reactome（up）.csv")



####下调DEG####################
#下调DEG
rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)

DOWN.DEG_rich <- all_DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  select(X) %>%
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
write.csv(kegg@result, "6.1 时序EC的KEGG（down）.csv")



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
write.csv(goMF@result, "6.2 时序EC的GO_MF（down）.csv")



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
write.csv(goCC@result, "6.3 时序EC的GO_CC（down）.csv")



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
write.csv(goBP@result, "6.4 时序EC的GO_BP（down）.csv")



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
write.csv(Reactome@result, "6.5 时序EC的Reactome（down）.csv")
