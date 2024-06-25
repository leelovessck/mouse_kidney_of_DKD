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
  if(!require(limma))install.packages("limma")
  if(!require(BiocParallel))install.packages("BiocParallel")
  if(!require(ReactomePA))BiocManager::install("ReactomePA")
  if(!require(DOSE))install.packages("DOSE")
}



####修改下载协议#####################
#修改下载协议
R.utils::setOption("clusterProfiler.download.method","auto")



####载入数据###############
#载入数据
ALL.DEG <- read.csv("./result/5.1.1 EC亚群DEG/全部EC的DEG.csv")

rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"
colnames(ALL.DEG)[3] <- "LogFC"
ALL.DEG <- ALL.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  select(SYMBOL, LogFC)
ALL.DEG_rich <- ALL.DEG$SYMBOL
entrezID <- bitr(ALL.DEG_rich,fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
ALL.DEG <- inner_join(ALL.DEG, entrezID,
                      by = "SYMBOL")
ALL.DEG<-ALL.DEG[,-1]
ALL.DEG<-na.omit(ALL.DEG)
ALL.DEG <- ALL.DEG %>%   
  arrange(desc(LogFC))

genelist = ALL.DEG[["LogFC"]]
names(genelist) = as.character(ALL.DEG[["ENTREZID"]])



####KEGG-GSEA###############
#KEGG-GSEA
KEGG_GSEA <- gseKEGG(geneList = genelist,
                     organism = "mmu",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     eps = 1e-10)
KEGG_GSEA_result <- KEGG_GSEA@result
write.csv(KEGG_GSEA_result, "1.全部EC的GSEA-KEGG.csv")

dotplot(KEGG_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到1.全部EC的GSEA-KEGG（点图）

significant_indices <- which(KEGG_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- KEGG_GSEA@result$Description[i]  
  plot1 <- gseaplot(KEGG_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.全部EC的GSEA-KEGG（各折线图）-无



####GO_GSEA###############
#GO_GSEA
GO_GSEA <- gseGO(genelist,
                 ont = "ALL",
                 OrgDb = 'org.Mm.eg.db',
                 keyType = "ENTREZID",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
GO_GSEA_result <- GO_GSEA@result
write.csv(GO_GSEA_result, "2.全部EC的GSEA-GO.csv")

dotplot(GO_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到2.全部EC的GSEA-GO（点图）

significant_indices <- which(GO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- GO_GSEA@result$Description[i]  
  plot1 <- gseaplot(GO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到2.全部EC的GSEA-GO（各折线图）


Reactome_GSEA <- gsePathway(genelist,
                            organism = "mouse",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500,
                            eps = 1e-10,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)
Reactome_GSEA_result <- Reactome_GSEA@result
write.csv(Reactome_GSEA_result, "3.全部EC的GSEA-Reactome.csv")

dotplot(Reactome_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到3.全部EC的GSEA-Reactome（点图）

significant_indices <- which(Reactome_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- Reactome_GSEA@result$Description[i]  
  plot1 <- gseaplot(Reactome_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到3.全部EC的GSEA-Reactome（各折线图）



####hallmark##################
#hallmark
hallmark <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_GSEA <- GSEA(genelist,
               TERM2GENE = hallmark,
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 1,
               pAdjustMethod = "BH")
H_GSEA_result <- H_GSEA@result
write.csv(H_GSEA_result, "4.全部EC的GSEA-hallmark.csv")

dotplot(H_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到4.全部EC的GSEA-hallmark（点图）

significant_indices <- which(H_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- H_GSEA@result$Description[i]  
  plot1 <- gseaplot(H_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.全部EC的GSEA-hallmark（各折线图）
