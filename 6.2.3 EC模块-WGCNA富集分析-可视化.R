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
all_kegg <- read.csv("./result/6.2.1 EC模块-WGCNA富集分析/1.EC模块的KEGG.csv")
all_gomf <- read.csv("./result/6.2.1 EC模块-WGCNA富集分析/2.EC模块的GO_MF.csv")
all_gocc <- read.csv("./result/6.2.1 EC模块-WGCNA富集分析/3.EC模块的GO_CC.csv")
all_gobp <- read.csv("./result/6.2.1 EC模块-WGCNA富集分析/4.EC模块的GO_BP.csv")
all_reactome <- read.csv("./result/6.2.1 EC模块-WGCNA富集分析/5.EC模块的Reactome.csv")



####可视化####################
#可视化
all_kegg_p <- all_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_kegg_p)
  ##得到1.EC模块_kegg

all_gomf_p <- all_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_gomf_p)
  ##得到2.EC模块_GOMF

all_gocc_p <- all_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_gocc_p)
  ##得到3.EC模块_GOCC

all_gobp_p <- all_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(all_gobp_p)
  ##得到4.EC模块_GOBP

all_reactome_p <- all_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_reactome_p)
  ##得到5.EC模块_reactome
