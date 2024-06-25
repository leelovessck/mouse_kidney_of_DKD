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
all_kegg <- read.csv("./result/6.1.1 EC-DEG富集分析/1.1 全部EC的KEGG（all）.csv")
all_gomf <- read.csv("./result/6.1.1 EC-DEG富集分析/1.2 全部EC的GO_MF（all）.csv")
all_gocc <- read.csv("./result/6.1.1 EC-DEG富集分析/1.3 全部EC的GO_CC（all）.csv")
all_gobp <- read.csv("./result/6.1.1 EC-DEG富集分析/1.4 全部EC的GO_BP（all）.csv")
all_reactome <- read.csv("./result/6.1.1 EC-DEG富集分析/1.5 全部EC的Reactome（all）.csv")

up_kegg <- read.csv("./result/6.1.1 EC-DEG富集分析/2.1 全部EC的KEGG（up）.csv")
up_gomf <- read.csv("./result/6.1.1 EC-DEG富集分析/2.2 全部EC的GO_MF（up）.csv")
up_gocc <- read.csv("./result/6.1.1 EC-DEG富集分析/2.3 全部EC的GO_CC（up）.csv")
up_gobp <- read.csv("./result/6.1.1 EC-DEG富集分析/2.4 全部EC的GO_BP（up）.csv")
up_reactome <- read.csv("./result/6.1.1 EC-DEG富集分析/2.5 全部EC的Reactome（up）.csv")

down_kegg <- read.csv("./result/6.1.1 EC-DEG富集分析/3.1 全部EC的KEGG（down）.csv")
down_gomf <- read.csv("./result/6.1.1 EC-DEG富集分析/3.2 全部EC的GO_MF（down）.csv")
down_gocc <- read.csv("./result/6.1.1 EC-DEG富集分析/3.3 全部EC的GO_CC（down）.csv")
down_gobp <- read.csv("./result/6.1.1 EC-DEG富集分析/3.4 全部EC的GO_BP（down）.csv")
down_reactome <- read.csv("./result/6.1.1 EC-DEG富集分析/3.5 全部EC的Reactome（down）.csv")



####可视化ALL####################
#可视化ALL
all_kegg_p <- all_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_kegg_p)
  ##得到1.1 all_kegg

all_gomf_p <- all_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_gomf_p)
  ##得到1.2 all_GOMF

all_gocc_p <- all_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_gocc_p)
  ##得到1.3 all_GOCC

all_gobp_p <- all_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(all_gobp_p)
  ##得到1.4 all_GOBP

all_reactome_p <- all_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(all_reactome_p)
  ##得到1.5 all_reactome



####可视化UP####################
#可视化UP
up_kegg_p <- up_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(up_kegg_p)
  ##得到2.1 up_kegg

up_gomf_p <- up_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(up_gomf_p)
  ##得到2.2 up_GOMF

up_gocc_p <- up_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(up_gocc_p)
  ##得到2.3 up_GOCC

up_gobp_p <- up_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(up_gobp_p)
  ##得到2.4 up_GOBP

up_reactome_p <- up_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(up_reactome_p)
  ##得到2.5 up_reactome



####可视化DOWN####################
#可视化DOWN
down_kegg_p <- down_kegg %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(down_kegg_p)
  ##得到3.1 down_kegg

down_gomf_p <- down_gomf %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(down_gomf_p)
  ##得到3.2 down_GOMF

down_gocc_p <- down_gocc %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(down_gocc_p)
  ##得到3.3 down_GOCC

down_gobp_p <- down_gobp %>% 
  filter(p.adjust < 0.05)
enrichmentNetwork(down_gobp_p)
  ##得到3.4 down_GOBP

down_reactome_p <- down_reactome %>% 
  filter(pvalue < 0.05)
enrichmentNetwork(down_reactome_p)
  ##得到3.5 down_reactome
