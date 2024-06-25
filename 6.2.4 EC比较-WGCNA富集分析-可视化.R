rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()

#打开必要的package
{
  if(!require(devtools))BiocManager::install("devtools")
  if(!require(aPEAR))install.packages("aPEAR")
  if(!require(clusterProfiler))install.packages("clusterProfiler")
  if(!require(org.Mm.eg.db))install.packages("org.Mm.eg.db")
  if(!require(DOSE))install.packages("DOSE")
  if(!require(ggplot2))install.packages("ggplot2")
}

####载入数据###############
#载入数据
dkd_kegg <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/1.1 EC-DKD的KEGG.csv")
dkd_gomf <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/1.2 EC-DKD的GO_MF.csv")
dkd_gocc <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/1.3 EC-DKD的GO_CC.csv")
dkd_gobp <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/1.4 EC-DKD的GO_BP.csv")
dkd_reactome <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/1.5 EC-DKD的Reactome.csv")

con_kegg <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/2.1 EC-CON的KEGG.csv")
con_gomf <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/2.2 EC-CON的GO_MF.csv")
con_gocc <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/2.3 EC-CON的GO_CC.csv")
con_gobp <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/2.4 EC-CON的GO_BP.csv")
con_reactome <- read.csv("./result/6.2.2 EC比较-WGCNA富集分析/2.5 EC-CON的Reactome.csv")



####可视化dkd####################
#可视化dkd
dkd_kegg_p <- dkd_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(dkd_kegg_p)
  ##得到1.1 dkd_kegg

dkd_gomf_p <- dkd_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(dkd_gomf_p)
  ##得到1.2 dkd_GOMF

dkd_gocc_p <- dkd_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(dkd_gocc_p)
  ##得到1.3 dkd_GOCC

dkd_gobp_p <- dkd_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(dkd_gobp_p)
  ##得到1.4 dkd_GOBP

dkd_reactome_p <- dkd_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(dkd_reactome_p)
  ##得到1.5 dkd_reactome



####可视化con####################
#可视化con
con_kegg_p <- con_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(con_kegg_p)
  ##得到2.1 con_kegg

con_gomf_p <- con_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(con_gomf_p)
  ##得到2.2 con_GOMF

con_gocc_p <- con_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(con_gocc_p)
  ##得到2.3 con_GOCC

con_gobp_p <- con_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(con_gobp_p)
  ##得到2.4 con_GOBP

con_reactome_p <- con_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(con_reactome_p)
  ##得到2.5 con_reactome
