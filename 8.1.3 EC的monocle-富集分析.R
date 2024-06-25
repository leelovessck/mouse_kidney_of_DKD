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
ALL.DEG <- read.csv("./result/8.2 EC的monocle-可视化/时序模块基因.csv")
rownames(ALL.DEG) <- ALL.DEG[,2]



####C1的时序基因####################
#C1的时序基因
C1.DEG_rich <- ALL.DEG %>%     
  filter(cluster == "cluster 1") %>%
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

C1.DEG_df <- bitr(C1.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
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
write.csv(kegg@result, "1.5 时序基因的Reactome（C1）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)



####C2的时序基因####################
#C2的时序基因
C2.DEG_rich <- ALL.DEG %>%     
  filter(cluster == "cluster 2") %>%
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

C2.DEG_df <- bitr(C2.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
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
write.csv(kegg@result, "2.5 时序基因的Reactome（C2）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)



####C3的时序基因####################
#C3的时序基因
C3.DEG_rich <- ALL.DEG %>%     
  filter(cluster == "cluster 3") %>%
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

C3.DEG_df <- bitr(C3.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
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
write.csv(kegg@result, "3.5 时序基因的Reactome（C3）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)



####C4的时序基因####################
#C4的时序基因
C4.DEG_rich <- ALL.DEG %>%     
  filter(cluster == "cluster 4") %>%
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

C4.DEG_df <- bitr(C4.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
C4.DEG_df <- C4.DEG_df %>% 
  distinct(SYMBOL, .keep_all = T)



####KEGG富集分析-C4#################
#KEGG富集分析-C4
kegg <- enrichKEGG(unique(C4.DEG_df$ENTREZID), 
                   organism ='mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = 0.2,
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = FALSE)
kegg <- setReadable(kegg,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "4.1 时序基因的KEGG（C4）.csv")



####GO（molecular function）-C4#################
#GO（molecular function）-C4
goMF <- enrichGO(C4.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goMF@result, "4.2 时序基因的GO_MF（C4）.csv")



####GO（cell component）-C4######################
#GO（cell component）-C4
goCC <- enrichGO(C4.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goCC@result, "4.3 时序基因的GO_CC（C4）.csv")



####GO（biological process）-C4##########################
#GO（biological process）-C4
goBP <- enrichGO(C4.DEG_df$SYMBOL, 
                 org.Mm.eg.db, 
                 keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 minGSSize = 10,
                 maxGSSize = 500, 
                 readable = FALSE, 
                 pool = FALSE)
write.csv(goBP@result, "4.4 时序基因的GO_BP（C4）.csv")



####Reactome-C4##########################
#Reactome-C4
Reactome <- enrichPathway(unique(C4.DEG_df$ENTREZID),
                          organism = "mouse",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Mm.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "4.5 时序基因的Reactome（C4）.csv")

rm(goBP, goCC, goMF, kegg, kegg_category, Reactome)
