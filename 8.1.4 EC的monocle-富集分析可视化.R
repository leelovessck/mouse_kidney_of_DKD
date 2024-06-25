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
C1_kegg <- read.csv("./result/8.3 EC的monocle-富集分析/1.1 时序基因的KEGG（C1）.csv")
C1_gomf <- read.csv("./result/8.3 EC的monocle-富集分析/1.2 时序基因的GO_MF（C1）.csv")
C1_gocc <- read.csv("./result/8.3 EC的monocle-富集分析/1.3 时序基因的GO_CC（C1）.csv")
C1_gobp <- read.csv("./result/8.3 EC的monocle-富集分析/1.4 时序基因的GO_BP（C1）.csv")
C1_reactome <- read.csv("./result/8.3 EC的monocle-富集分析/1.5 时序基因的Reactome（C1）.csv")

C2_kegg <- read.csv("./result/8.3 EC的monocle-富集分析/2.1 时序基因的KEGG（C2）.csv")
C2_gomf <- read.csv("./result/8.3 EC的monocle-富集分析/2.2 时序基因的GO_MF（C2）.csv")
C2_gocc <- read.csv("./result/8.3 EC的monocle-富集分析/2.3 时序基因的GO_CC（C2）.csv")
C2_gobp <- read.csv("./result/8.3 EC的monocle-富集分析/2.4 时序基因的GO_BP（C2）.csv")
C2_reactome <- read.csv("./result/8.3 EC的monocle-富集分析/2.5 时序基因的Reactome（C2）.csv")

C3_kegg <- read.csv("./result/8.3 EC的monocle-富集分析/3.1 时序基因的KEGG（C3）.csv")
C3_gomf <- read.csv("./result/8.3 EC的monocle-富集分析/3.2 时序基因的GO_MF（C3）.csv")
C3_gocc <- read.csv("./result/8.3 EC的monocle-富集分析/3.3 时序基因的GO_CC（C3）.csv")
C3_gobp <- read.csv("./result/8.3 EC的monocle-富集分析/3.4 时序基因的GO_BP（C3）.csv")
C3_reactome <- read.csv("./result/8.3 EC的monocle-富集分析/3.5 时序基因的Reactome（C3）.csv")

C4_kegg <- read.csv("./result/8.3 EC的monocle-富集分析/4.1 时序基因的KEGG（C4）.csv")
C4_gomf <- read.csv("./result/8.3 EC的monocle-富集分析/4.2 时序基因的GO_MF（C4）.csv")
C4_gocc <- read.csv("./result/8.3 EC的monocle-富集分析/4.3 时序基因的GO_CC（C4）.csv")
C4_gobp <- read.csv("./result/8.3 EC的monocle-富集分析/4.4 时序基因的GO_BP（C4）.csv")
C4_reactome <- read.csv("./result/8.3 EC的monocle-富集分析/4.5 时序基因的Reactome（C4）.csv")



####可视化C1####################
#可视化C1
C1_kegg_p <- C1_kegg %>% 
  filter(pvalue < 0.05)
C1p1 <- enrichmentNetwork(C1_kegg_p)
filename <- "1.1 C1_kegg.png"  
png(filename, width = 1000, height = 750)  
print(C1p1)  
dev.off()
  ##得到1.1 C1_kegg

C1_gomf_p <- C1_gomf %>% 
  filter(pvalue < 0.05)
C1p2 <- enrichmentNetwork(C1_gomf_p)
filename <- "1.2 C1_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(C1p2)  
dev.off()
  ##得到1.2 C1_GOMF

C1_gocc_p <- C1_gocc %>% 
  filter(pvalue < 0.05)
C1p3 <- enrichmentNetwork(C1_gocc_p)
filename <- "1.3 C1_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(C1p3)  
dev.off()
  ##得到1.3 C1_GOCC

C1_gobp_p <- C1_gobp %>% 
  filter(p.adjust < 0.05)
C1p4 <- enrichmentNetwork(C1_gobp_p)
filename <- "1.4 C1_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(C1p4)  
dev.off()
  ##得到1.4 C1_GOBP

C1_reactome_p <- C1_reactome %>% 
  filter(pvalue < 0.05)
C1p5 <- enrichmentNetwork(C1_reactome_p)
filename <- "1.5 C1_reactome.png"  
png(filename, width = 1000, height = 750)  
print(C1p5)  
dev.off()
  ##得到1.5 C1_reactome



####可视化C2####################
#可视化C2
C2_kegg_p <- C2_kegg %>% 
  filter(pvalue < 0.05)
C2p1 <- enrichmentNetwork(C2_kegg_p)
filename <- "2.1 C2_kegg.png"  
png(filename, width = 1000, height = 750)  
print(C2p1)  
dev.off()
##得到2.1 C2_kegg

C2_gomf_p <- C2_gomf %>% 
  filter(pvalue < 0.05)
C2p2 <- enrichmentNetwork(C2_gomf_p)
filename <- "2.2 C2_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(C2p2)  
dev.off()
##得到2.2 C2_GOMF

C2_gocc_p <- C2_gocc %>% 
  filter(pvalue < 0.05)
C2p3 <- enrichmentNetwork(C2_gocc_p)
filename <- "2.3 C2_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(C2p3)  
dev.off()
##得到2.3 C2_GOCC

C2_gobp_p <- C2_gobp %>% 
  filter(p.adjust < 0.05)
C2p4 <- enrichmentNetwork(C2_gobp_p)
filename <- "2.4 C2_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(C2p4)  
dev.off()
##得到2.4 C2_GOBP

C2_reactome_p <- C2_reactome %>% 
  filter(pvalue < 0.05)
C2p5 <- enrichmentNetwork(C2_reactome_p)
filename <- "2.5 C2_reactome.png"  
png(filename, width = 1000, height = 750)  
print(C2p5)  
dev.off()
##得到2.5 C2_reactome



####可视化C3####################
#可视化C3
C3_kegg_p <- C3_kegg %>% 
  filter(pvalue < 0.05)
C3p1 <- enrichmentNetwork(C3_kegg_p)
filename <- "3.1 C3_kegg.png"  
png(filename, width = 1000, height = 750)  
print(C3p1)  
dev.off()
##得到3.1 C3_kegg

C3_gomf_p <- C3_gomf %>% 
  filter(pvalue < 0.05)
C3p2 <- enrichmentNetwork(C3_gomf_p)
filename <- "3.2 C3_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(C3p2)  
dev.off()
##得到3.2 C3_GOMF

C3_gocc_p <- C3_gocc %>% 
  filter(pvalue < 0.05)
C3p3 <- enrichmentNetwork(C3_gocc_p)
filename <- "3.3 C3_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(C3p3)  
dev.off()
##得到3.3 C3_GOCC

C3_gobp_p <- C3_gobp %>% 
  filter(p.adjust < 0.05)
C3p4 <- enrichmentNetwork(C3_gobp_p)
filename <- "3.4 C3_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(C3p4)  
dev.off()
##得到3.4 C3_GOBP

C3_reactome_p <- C3_reactome %>% 
  filter(pvalue < 0.05)
C3p5 <- enrichmentNetwork(C3_reactome_p)
filename <- "3.5 C3_reactome.png"  
png(filename, width = 1000, height = 750)  
print(C3p5)  
dev.off()
##得到3.5 C3_reactome



####可视化C4####################
#可视化C4
C4_kegg_p <- C4_kegg %>% 
  filter(pvalue < 0.05)
C4p1 <- enrichmentNetwork(C4_kegg_p)
filename <- "4.1 C4_kegg.png"  
png(filename, width = 1000, height = 750)  
print(C4p1)  
dev.off()
##得到4.1 C4_kegg

C4_gomf_p <- C4_gomf %>% 
  filter(pvalue < 0.05)
C4p2 <- enrichmentNetwork(C4_gomf_p)
filename <- "4.2 C4_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(C4p2)  
dev.off()
##得到4.2 C4_GOMF

C4_gocc_p <- C4_gocc %>% 
  filter(pvalue < 0.05)
C4p3 <- enrichmentNetwork(C4_gocc_p)
filename <- "4.3 C4_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(C4p3)  
dev.off()
##得到4.3 C4_GOCC

C4_gobp_p <- C4_gobp %>% 
  filter(p.adjust < 0.05)
C4p4 <- enrichmentNetwork(C4_gobp_p)
filename <- "4.4 C4_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(C4p4)  
dev.off()
##得到4.4 C4_GOBP

C4_reactome_p <- C4_reactome %>% 
  filter(pvalue < 0.05)
C4p5 <- enrichmentNetwork(C4_reactome_p)
filename <- "4.5 C4_reactome.png"  
png(filename, width = 1000, height = 750)  
print(C4p5)  
dev.off()
##得到4.5 C4_reactome


