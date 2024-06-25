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
C1_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/1.1 时序基因的KEGG（C1）.csv")
C1_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/1.2 时序基因的GO_MF（C1）.csv")
C1_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/1.3 时序基因的GO_CC（C1）.csv")
C1_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/1.4 时序基因的GO_BP（C1）.csv")
C1_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/1.5 时序基因的Reactome（C1）.csv")

C2_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/2.1 时序基因的KEGG（C2）.csv")
C2_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/2.2 时序基因的GO_MF（C2）.csv")
C2_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/2.3 时序基因的GO_CC（C2）.csv")
C2_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/2.4 时序基因的GO_BP（C2）.csv")
C2_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/2.5 时序基因的Reactome（C2）.csv")

C3_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/3.1 时序基因的KEGG（C3）.csv")
C3_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/3.2 时序基因的GO_MF（C3）.csv")
C3_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/3.3 时序基因的GO_CC（C3）.csv")
C3_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/3.4 时序基因的GO_BP（C3）.csv")
C3_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/3.5 时序基因的Reactome（C3）.csv")

ALL_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/4.1 时序EC的KEGG（all）.csv")
ALL_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/4.2 时序EC的GO_MF（all）.csv")
ALL_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/4.3 时序EC的GO_CC（all）.csv")
ALL_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/4.4 时序EC的GO_BP（all）.csv")
ALL_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/4.5 时序EC的Reactome（all）.csv")

UP_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/5.1 时序EC的KEGG（up）.csv")
UP_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/5.2 时序EC的GO_MF（up）.csv")
UP_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/5.3 时序EC的GO_CC（up）.csv")
UP_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/5.4 时序EC的GO_BP（up）.csv")
UP_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/5.5 时序EC的Reactome（up）.csv")

DOWN_kegg <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/6.1 时序EC的KEGG（down）.csv")
DOWN_gomf <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/6.2 时序EC的GO_MF（down）.csv")
DOWN_gocc <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/6.3 时序EC的GO_CC（down）.csv")
DOWN_gobp <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/6.4 时序EC的GO_BP（down）.csv")
DOWN_reactome <- read.csv("./result/8.5.2 DKD EC的阶段-富集分析/6.5 时序EC的Reactome（down）.csv")



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



####可视化all DEG####################
#可视化all DEG
ALL_kegg_p <- ALL_kegg %>% 
  filter(pvalue < 0.05)
ALLp1 <- enrichmentNetwork(ALL_kegg_p)
filename <- "4.1 ALL_kegg.png"  
png(filename, width = 1000, height = 750)  
print(ALLp1)  
dev.off()
##得到4.1 ALL_kegg

ALL_gomf_p <- ALL_gomf %>% 
  filter(pvalue < 0.05)
ALLp2 <- enrichmentNetwork(ALL_gomf_p)
filename <- "4.2 ALL_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(ALLp2)  
dev.off()
##得到4.2 ALL_GOMF

ALL_gocc_p <- ALL_gocc %>% 
  filter(pvalue < 0.05)
ALLp3 <- enrichmentNetwork(ALL_gocc_p)
filename <- "4.3 ALL_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(ALLp3)  
dev.off()
##得到4.3 ALL_GOCC

ALL_gobp_p <- ALL_gobp %>% 
  filter(p.adjust < 0.05)
ALLp4 <- enrichmentNetwork(ALL_gobp_p)
filename <- "4.4 ALL_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(ALLp4)  
dev.off()
##得到4.4 ALL_GOBP

ALL_reactome_p <- ALL_reactome %>% 
  filter(pvalue < 0.05)
ALLp5 <- enrichmentNetwork(ALL_reactome_p)
filename <- "4.5 ALL_reactome.png"  
png(filename, width = 1000, height = 750)  
print(ALLp5)  
dev.off()
##得到4.5 ALL_reactome



####可视化UP DEG####################
#可视化UP DEG
UP_kegg_p <- UP_kegg %>% 
  filter(pvalue < 0.05)
UPp1 <- enrichmentNetwork(UP_kegg_p)
filename <- "5.1 UP_kegg.png"  
png(filename, width = 1000, height = 750)  
print(UPp1)  
dev.off()
##得到5.1 UP_kegg

UP_gomf_p <- UP_gomf %>% 
  filter(pvalue < 0.05)
UPp2 <- enrichmentNetwork(UP_gomf_p)
filename <- "5.2 UP_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(UPp2)  
dev.off()
##得到5.2 UP_GOMF

UP_gocc_p <- UP_gocc %>% 
  filter(pvalue < 0.05)
UPp3 <- enrichmentNetwork(UP_gocc_p)
filename <- "5.3 UP_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(UPp3)  
dev.off()
##得到5.3 UP_GOCC

UP_gobp_p <- UP_gobp %>% 
  filter(p.adjust < 0.05)
UPp4 <- enrichmentNetwork(UP_gobp_p)
filename <- "5.4 UP_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(UPp4)  
dev.off()
##得到5.4 UP_GOBP

UP_reactome_p <- UP_reactome %>% 
  filter(pvalue < 0.05)
UPp5 <- enrichmentNetwork(UP_reactome_p)
filename <- "5.5 UP_reactome.png"  
png(filename, width = 1000, height = 750)  
print(UPp5)  
dev.off()
##得到5.5 UP_reactome



####可视化DOWN DEG####################
#可视化DOWN DEG
DOWN_kegg_p <- DOWN_kegg %>% 
  filter(pvalue < 0.05)
DOWNp1 <- enrichmentNetwork(DOWN_kegg_p)
filename <- "6.1 DOWN_kegg.png"  
png(filename, width = 1000, height = 750)  
print(DOWNp1)  
dev.off()
##得到6.1 DOWN_kegg

DOWN_gomf_p <- DOWN_gomf %>% 
  filter(pvalue < 0.05)
DOWNp2 <- enrichmentNetwork(DOWN_gomf_p)
filename <- "6.2 DOWN_GOMF.png"  
png(filename, width = 1000, height = 750)  
print(DOWNp2)  
dev.off()
##得到6.2 DOWN_GOMF

DOWN_gocc_p <- DOWN_gocc %>% 
  filter(pvalue < 0.05)
DOWNp3 <- enrichmentNetwork(DOWN_gocc_p)
filename <- "6.3 DOWN_GOCC.png"  
png(filename, width = 1000, height = 750)  
print(DOWNp3)  
dev.off()
##得到6.3 DOWN_GOCC

DOWN_gobp_p <- DOWN_gobp %>% 
  filter(p.adjust < 0.05)
DOWNp4 <- enrichmentNetwork(DOWN_gobp_p)
filename <- "6.4 DOWN_GOBP.png"  
png(filename, width = 1000, height = 750)  
print(DOWNp4)  
dev.off()
##得到6.4 DOWN_GOBP

DOWN_reactome_p <- DOWN_reactome %>% 
  filter(pvalue < 0.05)
DOWNp5 <- enrichmentNetwork(DOWN_reactome_p)
filename <- "6.5 DOWN_reactome.png"  
png(filename, width = 1000, height = 750)  
print(DOWNp5)  
dev.off()
##得到6.5 DOWN_reactome
