rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
  if(!require(ggsignif))install.packages("ggsignif")
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
klk_rows <- grepl("Klk", expression_CONEC$gene_id, ignore.case = FALSE)  
KLK_CONEC <- expression_CONEC[klk_rows, ]  
non_zero_counts <- rowSums(KLK_CONEC[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_CONEC <- KLK_CONEC[rows_to_keep, ] 
write.csv(KLK_CONEC, "KLK_CONEC.csv")

express_DKDEC <- GetAssayData(DKD_EC, layer = "data")
expression_DKDEC <- as.data.frame(express_DKDEC)
expression_DKDEC <- expression_DKDEC %>% 
  rownames_to_column("gene_id")
klk_rows <- grepl("Klk", expression_DKDEC$gene_id, ignore.case = FALSE)  
KLK_DKDEC <- expression_DKDEC[klk_rows, ]  
non_zero_counts <- rowSums(KLK_DKDEC[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_DKDEC <- KLK_DKDEC[rows_to_keep, ] 
write.csv(KLK_DKDEC, "KLK_DKDEC.csv")



####读取长数据基因（EC）###############
#读取长数据基因（EC）
mouse <- read.csv("./result/9.1.1 KLK-基因/长-KLK_EC.csv")

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
filename <- "1.1 EC的Klk表达量总览.png"  
png(filename, width = 1000, height = 1000)  
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
filename <- "1.2 EC的Klk表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 1000)  
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
klk_rows <- grepl("Klk", expression_CONPOD$gene_id, ignore.case = FALSE)  
KLK_CONPOD <- expression_CONPOD[klk_rows, ]  
non_zero_counts <- rowSums(KLK_CONPOD[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_CONPOD <- KLK_CONPOD[rows_to_keep, ] 
write.csv(KLK_CONPOD, "KLK_CONPOD.csv")

express_DKDPOD <- GetAssayData(DKD_POD, layer = "data")
expression_DKDPOD <- as.data.frame(express_DKDPOD)
expression_DKDPOD <- expression_DKDPOD %>% 
  rownames_to_column("gene_id")
klk_rows <- grepl("Klk", expression_DKDPOD$gene_id, ignore.case = FALSE)  
KLK_DKDPOD <- expression_DKDPOD[klk_rows, ]  
non_zero_counts <- rowSums(KLK_DKDPOD[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_DKDPOD <- KLK_DKDPOD[rows_to_keep, ] 
write.csv(KLK_DKDPOD, "KLK_DKDPOD.csv")



####读取长数据基因（POD）###############
#读取长数据基因（POD）
mouse <- read.csv("./result/9.1.1 KLK-基因/长-KLK_POD.csv")

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
filename <- "2.1 POD的Klk表达量总览.png"  
png(filename, width = 1000, height = 1000)  
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
filename <- "2.2 POD的Klk表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 1000)  
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
klk_rows <- grepl("Klk", expression_CONMC$gene_id, ignore.case = FALSE)  
KLK_CONMC <- expression_CONMC[klk_rows, ]  
non_zero_counts <- rowSums(KLK_CONMC[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_CONMC <- KLK_CONMC[rows_to_keep, ] 
write.csv(KLK_CONMC, "KLK_CONMC.csv")

express_DKDMC <- GetAssayData(DKD_MC, layer = "data")
expression_DKDMC <- as.data.frame(express_DKDMC)
expression_DKDMC <- expression_DKDMC %>% 
  rownames_to_column("gene_id")
klk_rows <- grepl("Klk", expression_DKDMC$gene_id, ignore.case = FALSE)  
KLK_DKDMC <- expression_DKDMC[klk_rows, ]  
non_zero_counts <- rowSums(KLK_DKDMC[,-1] != 0)  
rows_to_keep <- non_zero_counts > 0  
KLK_DKDMC <- KLK_DKDMC[rows_to_keep, ] 
write.csv(KLK_DKDMC, "KLK_DKDMC.csv")



####读取长数据基因（MC）###############
#读取长数据基因（MC）
mouse <- read.csv("./result/9.1.1 KLK-基因/长-KLK_MC.csv")

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
filename <- "3.1 MC的Klk表达量总览.png"  
png(filename, width = 1000, height = 1000)  
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
filename <- "3.2 MC的Klk表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 1000)  
print(p1)  
dev.off()
