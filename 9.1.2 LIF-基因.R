rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")



####提取基因（EC）############################
#提取基因（EC）
EC <- subset(mouse_all, usetype == "EC")
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
table(EC$usetype)

CON_EC <- subset(EC, group == "CON")
table(CON_EC$usetype)
DKD_EC <- subset(EC, group == "DKD")
table(DKD_EC$usetype)

express_CONEC <- GetAssayData(CON_EC, layer = "data")
expression_CONEC <- as.data.frame(express_CONEC)
expression_CONEC <- expression_CONEC %>% 
  rownames_to_column("gene_id")
LIF_CONEC <- expression_CONEC %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_CONEC, "LIF_CONEC.csv")

express_DKDEC <- GetAssayData(DKD_EC, layer = "data")
expression_DKDEC <- as.data.frame(express_DKDEC)
expression_DKDEC <- expression_DKDEC %>% 
  rownames_to_column("gene_id")
LIF_DKDEC <- expression_DKDEC %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_DKDEC, "LIF_DKDEC.csv")



####读取长数据基因（EC）###############
#读取长数据基因（EC）
mouse <- read.csv("./result/9.1.2 LIF-基因/长-LIF_EC.csv")

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "free")
p1
filename <- "1.1 EC的Lif表达量总览.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "fixed")
p1
filename <- "1.2 EC的Lif表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()



####提取基因（POD）############################
#提取基因（POD）
POD <- subset(mouse_all, usetype == "POD")
POD$usetype <- droplevels(POD$usetype,
                         exclude = setdiff(
                           levels(POD$usetype),
                           unique(POD$usetype)))
table(POD$usetype)

CON_POD <- subset(POD, group == "CON")
table(CON_POD$usetype)
DKD_POD <- subset(POD, group == "DKD")
table(DKD_POD$usetype)

express_CONPOD <- GetAssayData(CON_POD, layer = "data")
expression_CONPOD <- as.data.frame(express_CONPOD)
expression_CONPOD <- expression_CONPOD %>% 
  rownames_to_column("gene_id")
LIF_CONPOD <- expression_CONPOD %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_CONPOD, "LIF_CONPOD.csv")

express_DKDPOD <- GetAssayData(DKD_POD, layer = "data")
expression_DKDPOD <- as.data.frame(express_DKDPOD)
expression_DKDPOD <- expression_DKDPOD %>% 
  rownames_to_column("gene_id")
LIF_DKDPOD <- expression_DKDPOD %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_DKDPOD, "LIF_DKDPOD.csv")



####读取长数据基因（POD）###############
#读取长数据基因（POD）
mouse <- read.csv("./result/9.1.2 LIF-基因/长-LIF_POD.csv")

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "free")
p1
filename <- "2.1 POD的Lif表达量总览.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "fixed")
p1
filename <- "2.2 POD的Lif表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()



####提取基因（MC）############################
#提取基因（MC）
MC <- subset(mouse_all, usetype == "MC")
MC$usetype <- droplevels(MC$usetype,
                          exclude = setdiff(
                            levels(MC$usetype),
                            unique(MC$usetype)))
table(MC$usetype)

CON_MC <- subset(MC, group == "CON")
table(CON_MC$usetype)
DKD_MC <- subset(MC, group == "DKD")
table(DKD_MC$usetype)

express_CONMC <- GetAssayData(CON_MC, layer = "data")
expression_CONMC <- as.data.frame(express_CONMC)
expression_CONMC <- expression_CONMC %>% 
  rownames_to_column("gene_id")
LIF_CONMC <- expression_CONMC %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_CONMC, "LIF_CONMC.csv")

express_DKDMC <- GetAssayData(DKD_MC, layer = "data")
expression_DKDMC <- as.data.frame(express_DKDMC)
expression_DKDMC <- expression_DKDMC %>% 
  rownames_to_column("gene_id")
LIF_DKDMC <- expression_DKDMC %>% 
  filter(gene_id == "Lif" | gene_id == "Lifr")
write.csv(LIF_DKDMC, "LIF_DKDMC.csv")



####读取长数据基因（MC）###############
#读取长数据基因（MC）
mouse <- read.csv("./result/9.1.2 LIF-基因/长-LIF_MC.csv")

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "free")
p1
filename <- "3.1 MC的Lif表达量总览.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()

p1 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "fixed")
p1
filename <- "3.2 MC的Lif表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()
