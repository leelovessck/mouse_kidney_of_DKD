rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()

#打开必要的package
{
if(!require(NMF))install.packages('NMF')
if(!require(circlize))devtools::install_github("jokergoo/circlize")
if(!require(ComplexHeatmap))devtools::install_github("jokergoo/ComplexHeatmap")
if(!require(CellChat))devtools::install_github("jinworks/CellChat")
if(!require(Seurat))install.packages("Seurat")
if(!require(patchwork))install.packages("patchwork")
if(!require(cowplot))install.packages("cowplot")
if(!require(dplyr))install.packages("dplyr")
if(!require(ggalluvial))install.packages('ggalluvial')
if(!require(igraph))install.packages('igraph')
if(!require(forcats))install.packages("forcats")
if(!require(aplot))install.packages("aplot")
}



####数据准备############################
#数据准备
con.chat <- readRDS("./data source/cellchat/control（全部）.rds")
dkd.chat <- readRDS("./data source/cellchat/DKD（全部）.rds")
object.list <- list(CON = con.chat, DKD = dkd.chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
levels(object.list[[1]]@idents) 
levels(object.list[[2]]@idents) 
group.cellType <- levels(object.list[[1]]@idents)
group.cellType <- factor(group.cellType, levels = levels(object.list[[1]]@idents))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})#合并细胞类型
con.chat <- netAnalysis_computeCentrality(con.chat) 
dkd.chat <- netAnalysis_computeCentrality(dkd.chat) 
object.list <- list(CON = con.chat, DKD = dkd.chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



####差异计算#########################
#差异计算
pos.dataset = "DKD"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)



####提取DKD中改变的配受体对########################
#提取DKD中改变的配受体对
net <- netMappingDEG(cellchat, 
                     features.name = features.name)
ligand.up <- subsetCommunication(cellchat, 
                                 net = net, 
                                 datasets = "DKD",
                                 ligand.logFC = 0.2, 
                                 receptor.logFC = NULL)
receptor.up <- subsetCommunication(cellchat, 
                                   net = net, 
                                   datasets = "DKD",
                                   ligand.logFC = NULL, 
                                   receptor.logFC = 0.2)
ligand.down <- subsetCommunication(cellchat, 
                                   net = net, 
                                   datasets = "DKD",
                                   ligand.logFC = -0.2, 
                                   receptor.logFC = NULL)
receptor.down <- subsetCommunication(cellchat, 
                                     net = net, 
                                     datasets = "DKD",
                                     ligand.logFC = NULL, 
                                     receptor.logFC = -0.2)



####EC的上调配体气泡图（source）#####################
#EC的上调配体气泡图（source）
gc.s.lup <- ligand.up %>%  
  filter(source == "EC")

unique(gc.s.lup$pathway_name)
gc.s.lup.use = gc.s.lup[, "interaction_name", drop = F]

gc.s.lup$interaction_name <- as.factor(gc.s.lup$interaction_name)
gc.s.lup$interaction_name <- fct_inorder(gc.s.lup$interaction_name)
gc.s.lup$interaction_name_2 <- as.factor(gc.s.lup$interaction_name_2)
gc.s.lup$interaction_name_2 <- fct_inorder(gc.s.lup$interaction_name_2)
gc.s.lup$pathway_name <- as.factor(gc.s.lup$pathway_name)
gc.s.lup$pathway_name <- fct_inorder(gc.s.lup$pathway_name)
gc.s.lup <- gc.s.lup %>% 
  mutate(p="")
gc.s.lup$signway <- paste(gc.s.lup$source,
                          "  ->  ", 
                          gc.s.lup$target) 

table(unique(gc.s.lup$interaction_name_2))
summary_data <- gc.s.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Ligand of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "#FFB6C1",
                        high = "#FF0000",
                        midpoint = 0) +
  geom_hline(yintercept=c(2.5, 3.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("DHEAS" = 1.5, 
                     "ADGRA" = 3, 
                     "ADGRE" = 4.5) 
gene_group <- ggplot(gc.s.lup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group

p1%>%insert_left(gene_group, width = 0.3)
  ##得到1.1 EC的source上调配体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.lup.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC的source上调配体气泡图B



####EC的上调配体气泡图（target）#####################
#EC的上调配体气泡图（target）
gc.t.lup <- ligand.up %>%  
  filter(target == "EC")

gc.t.lup <- gc.t.lup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.lup$pathway_name)
gc.t.lup.use = gc.t.lup[, "interaction_name", drop = F]
gc.t.lup$interaction_name <- as.factor(gc.t.lup$interaction_name)
gc.t.lup$interaction_name <- fct_inorder(gc.t.lup$interaction_name)
gc.t.lup$interaction_name_2 <- as.factor(gc.t.lup$interaction_name_2)
gc.t.lup$interaction_name_2 <- fct_inorder(gc.t.lup$interaction_name_2)
gc.t.lup$pathway_name <- as.factor(gc.t.lup$pathway_name)
gc.t.lup$pathway_name <- fct_inorder(gc.t.lup$pathway_name)
gc.t.lup <- gc.t.lup %>% 
  mutate(p="")

gc.t.lup$signway <- paste(gc.t.lup$source,
                          "  ->  ", 
                          gc.t.lup$target) 

table(unique(gc.t.lup$interaction_name_2))
summary_data <- gc.t.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.t.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Ligand of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low="#FFB6C1",
                        high="#FF0000") +
  geom_hline(yintercept=c(1.5, 3.5, 5.5, 7.5, 11.5, 12.5, 13.5, 15.5, 16.5, 
                          18.5, 19.5, 20.5, 21.5, 25.5, 26.5, 28.5, 30.5, 31.5, 
                          37.5, 38.5, 43.5, 46.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("AGRN" = 1, 
                     "ADGRE" = 2.5, 
                     "ADGRG" = 4.5, 
                     "CSF" = 6.5, 
                     "COLLAGEN" = 9.5, 
                     "CEACAM" = 12, 
                     "DHEAS" = 13, 
                     "EGF" = 14.5, 
                     "EPHB" = 16, 
                     "FGF" = 17.5, 
                     "GAP" = 19, 
                     "IGF" = 20, 
                     "JAM" = 21, 
                     "LAMININ" = 23.5, 
                     "MMP" = 26, 
                     "NECTIN" = 27.5, 
                     "NOTCH" = 29.5, 
                     "PECAM2" = 31, 
                     "SPP1" = 34.5, 
                     "TGFb" = 38, 
                     "THBS" = 41, 
                     "VEGF" = 45, 
                     "VCAM" = 47.5)  
gene_group <- ggplot(gc.t.lup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group
p1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.3 EC的target上调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.lup.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 3.5, 5.5, 7.5, 11.5, 12.5, 13.5, 15.5, 16.5, 
                          18.5, 19.5, 20.5, 21.5, 25.5, 26.5, 28.5, 30.5, 31.5, 
                          37.5, 38.5, 43.5, 46.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC的target上调配体气泡图B



####EC的上调受体气泡图（source）#####################
#EC的上调受体气泡图（source）
gc.s.rup <- receptor.up %>%  
  filter(source == "EC")

gc.s.rup <- gc.s.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.s.rup$pathway_name)
gc.s.rup.use = gc.s.rup[, "interaction_name", drop = F]
gc.s.rup$interaction_name <- as.factor(gc.s.rup$interaction_name)
gc.s.rup$interaction_name <- fct_inorder(gc.s.rup$interaction_name)
gc.s.rup$interaction_name_2 <- as.factor(gc.s.rup$interaction_name_2)
gc.s.rup$interaction_name_2 <- fct_inorder(gc.s.rup$interaction_name_2)
gc.s.rup$pathway_name <- as.factor(gc.s.rup$pathway_name)
gc.s.rup$pathway_name <- fct_inorder(gc.s.rup$pathway_name)
gc.s.rup <- gc.s.rup %>% 
  mutate(p="")
gc.s.rup$signway <- paste(gc.s.rup$source,
                          "  ->  ", 
                          gc.s.rup$target) 

table(unique(gc.s.rup$interaction_name_2))
summary_data <- gc.s.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.s.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Receptor of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "#FFB6C1",
                        high = "#FF0000") +
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 9.5, 10.5, 11.5, 13.5, 15.5, 16.5, 
                          17.5, 18.5, 19.5, 20.5, 21.5, 25.5, 26.5, 38.5, 41.5,
                          42.5, 44.5, 46.5, 47.5, 48.5, 50.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPTL" = 1, 
                     "AGRN" = 2, 
                     "ADGRG" = 3, 
                     "COLLAGEN" = 6.5, 
                     "CD34" = 10, 
                     "CEACAM" = 11, 
                     "EPHA" = 12.5, 
                     "EPHB" = 14.5, 
                     "FGF" = 16, 
                     "GAS" = 17, 
                     "GALECTIN" = 18, 
                     "GAP" = 19, 
                     "HSPG" = 20, 
                     "ICAM" = 21, 
                     "JAM" = 23.5, 
                     "KLK" = 26, 
                     "LAMININ" = 32.5, 
                     "NOTCH" = 40, 
                     "Netrin" = 42, 
                     "PDGF" = 43.5, 
                     "PTN" = 45.5, 
                     "PROS" = 47, 
                     "SEMA5" = 48, 
                     "SEMA6" = 49.5, 
                     "VCAM" = 51.5) 
gene_group <- ggplot(gc.s.rup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group

p1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.1 EC的source上调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rup.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 9.5, 10.5, 11.5, 13.5, 15.5, 16.5, 
                          17.5, 18.5, 19.5, 20.5, 21.5, 25.5, 26.5, 38.5, 41.5,
                          42.5, 44.5, 46.5, 47.5, 48.5, 50.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC的source上调受体气泡图B



####EC的上调受体气泡图（target）#####################
#EC的上调受体气泡图（target）
gc.t.rup <- receptor.up %>%  
  filter(target == "EC")

gc.t.rup <- gc.t.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.rup$pathway_name)
gc.t.rup.use = gc.t.rup[, "interaction_name", drop = F]
gc.t.rup$interaction_name <- as.factor(gc.t.rup$interaction_name)
gc.t.rup$interaction_name <- fct_inorder(gc.t.rup$interaction_name)
gc.t.rup$interaction_name_2 <- as.factor(gc.t.rup$interaction_name_2)
gc.t.rup$interaction_name_2 <- fct_inorder(gc.t.rup$interaction_name_2)
gc.t.rup$pathway_name <- as.factor(gc.t.rup$pathway_name)
gc.t.rup$pathway_name <- fct_inorder(gc.t.rup$pathway_name)
gc.t.rup <- gc.t.rup %>% 
  mutate(p="")

gc.t.rup$signway <- paste(gc.t.rup$source,
                          "  ->  ", 
                          gc.t.rup$target) 

table(unique(gc.t.rup$interaction_name_2))
summary_data <- gc.t.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.t.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Receptor of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='#FFB6C1',
                        high='#FF0000',
                        midpoint = 0.8) +
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPTL" = 1, 
                     "ADGRG" = 2.5, 
                     "FN1" = 4, 
                     "IGF" = 5, 
                     "MMP" = 6, 
                     "PTN" = 7, 
                     "SPP1" = 8, 
                     "TENASCIN" = 9, 
                     "THBS" = 11.5)  
gene_group <- ggplot(gc.t.rup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group
p1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.3 EC的target上调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rup.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC的target上调受体气泡图B



####EC的下调配体气泡图（source）#####################
#EC的下调配体气泡图（source）
gc.s.ldown <- ligand.down %>%  
  filter(source == "EC")

gc.s.ldown <- gc.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 
 
unique(gc.s.ldown$pathway_name)
gc.s.ldown.use = gc.s.ldown[, "interaction_name", drop = F]
gc.s.ldown$interaction_name <- as.factor(gc.s.ldown$interaction_name)
gc.s.ldown$interaction_name <- fct_inorder(gc.s.ldown$interaction_name)
gc.s.ldown$interaction_name_2 <- as.factor(gc.s.ldown$interaction_name_2)
gc.s.ldown$interaction_name_2 <- fct_inorder(gc.s.ldown$interaction_name_2)
gc.s.ldown$pathway_name <- as.factor(gc.s.ldown$pathway_name)
gc.s.ldown$pathway_name <- fct_inorder(gc.s.ldown$pathway_name)
gc.s.ldown <- gc.s.ldown %>% 
  mutate(p="")
gc.s.ldown$signway <- paste(gc.s.ldown$source,
                          "  ->  ", 
                          gc.s.ldown$target) 

table(unique(gc.s.ldown$interaction_name_2))
summary_data <- gc.s.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.s.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Ligand of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low="#000080",
                        high="#ADD8E6") +
  geom_hline(yintercept=c(2.5, 3.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("GAS" = 1.5, 
                     "KIT" = 3, 
                     "SEMA6" = 4.5) 
gene_grodown <- ggplot(gc.s.ldown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown

p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.1 EC-GC的source下调配体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.ldown.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-GC的source下调配体气泡图B



####EC的下调配体气泡图（target）#####################
#EC的下调配体气泡图（target）
gc.t.ldown <- ligand.down %>%  
  filter(target == "EC")

gc.t.ldown <- gc.t.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.ldown$pathway_name)
gc.t.ldown.use = gc.t.ldown[, "interaction_name", drop = F]
gc.t.ldown$interaction_name <- as.factor(gc.t.ldown$interaction_name)
gc.t.ldown$interaction_name <- fct_inorder(gc.t.ldown$interaction_name)
gc.t.ldown$interaction_name_2 <- as.factor(gc.t.ldown$interaction_name_2)
gc.t.ldown$interaction_name_2 <- fct_inorder(gc.t.ldown$interaction_name_2)
gc.t.ldown$pathway_name <- as.factor(gc.t.ldown$pathway_name)
gc.t.ldown$pathway_name <- fct_inorder(gc.t.ldown$pathway_name)
gc.t.ldown <- gc.t.ldown %>% 
  mutate(p="")

gc.t.ldown$signway <- paste(gc.t.ldown$source,
                          "  ->  ", 
                          gc.t.ldown$target) 

table(unique(gc.t.ldown$interaction_name_2))
summary_data <- gc.t.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.t.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Ligand of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(1.5, 2.5, 5.5, 14.5, 15.5, 16.5, 17.5, 19.5, 20.5, 
                          21.5, 23.5, 24.5, 25.5, 27.5, 28.5, 29.5, 30.5, 31.5,
                          32.5, 36.5, 37.5, 43.5, 46.5, 47.5, 54.5, 56.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPTL" = 1, 
                     "ANGPT" = 2, 
                     "ADGRE" = 4, 
                     "COLLAGEN" = 10, 
                     "CypA" = 15, 
                     "CD45" = 16, 
                     "CSPG4" = 17, 
                     "CEACAM" = 18.5, 
                     "EGF" = 20, 
                     "ESAM" = 21, 
                     "FGF" = 22.5, 
                     "GAS" = 24, 
                     "GRN" = 25, 
                     "GALECTIN" = 26.5, 
                     "IL2" = 28, 
                     "KIT" = 29, 
                     "LAMININ" = 30, 
                     "MMP" = 31, 
                     "PDGF" = 32, 
                     "PTN" = 34.5, 
                     "PLAU" = 37, 
                     "SPP1" = 40.5, 
                     "SEMA3" = 45, 
                     "SEMA6" = 47, 
                     "VEGF" = 51, 
                     "VCAM" = 55.5, 
                     "VISTA" = 57)  
gene_grodown <- ggplot(gc.t.ldown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown
p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.3 EC的target下调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.ldown.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 5.5, 14.5, 15.5, 16.5, 17.5, 19.5, 20.5, 
                          21.5, 23.5, 24.5, 25.5, 27.5, 28.5, 29.5, 30.5, 31.5,
                          32.5, 36.5, 37.5, 43.5, 46.5, 47.5, 54.5, 56.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC的target下调配体气泡图B



####EC的下调受体气泡图（source）#####################
#EC的下调受体气泡图（source）
gc.s.rdown <- receptor.down %>%  
  filter(source == "EC")

gc.s.rdown <- gc.s.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.s.rdown$pathway_name)
gc.s.rdown.use = gc.s.rdown[, "interaction_name", drop = F]
gc.s.rdown$interaction_name <- as.factor(gc.s.rdown$interaction_name)
gc.s.rdown$interaction_name <- fct_inorder(gc.s.rdown$interaction_name)
gc.s.rdown$interaction_name_2 <- as.factor(gc.s.rdown$interaction_name_2)
gc.s.rdown$interaction_name_2 <- fct_inorder(gc.s.rdown$interaction_name_2)
gc.s.rdown$pathway_name <- as.factor(gc.s.rdown$pathway_name)
gc.s.rdown$pathway_name <- fct_inorder(gc.s.rdown$pathway_name)
gc.s.rdown <- gc.s.rdown %>% 
  mutate(p="")
gc.s.rdown$signway <- paste(gc.s.rdown$source,
                          "  ->  ", 
                          gc.s.rdown$target) 

table(unique(gc.s.rdown$interaction_name_2))
summary_data <- gc.s.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.s.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Receptor of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5, 25.5, 27.5, 28.5,
                          29.5, 31.5, 33.5, 34.5, 37.5, 41.5, 42.5, 67.5, 68.5,
                          69.5, 71.5, 72.5, 77.5, 79.5, 81.5, 85.5, 86.5, 88.5,
                          89.5, 90.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, 
                     "AGRN" = 2, 
                     "APP" = 3, 
                     "ADGRE" = 4.5, 
                     "ADGRG" = 6, 
                     "CSF" = 7.5, 
                     "COMPLEMENT" = 9, 
                     "COLLAGEN" = 17.5, 
                     "CEACAM" = 26.5, 
                     "EPHB" = 28, 
                     "ESAM" = 29, 
                     "GAS" = 30.5, 
                     "GALECTIN" = 32.5, 
                     "HSPG" = 34, 
                     "ICAM" = 36, 
                     "JAM" = 39.5, 
                     "KLK" = 42, 
                     "LAMININ" = 55, 
                     "MIF" = 68, 
                     "MMP" = 69, 
                     "NOTCH" = 70.5, 
                     "Netrin" = 72, 
                     "PDGF" = 75, 
                     "PTN" = 78.5, 
                     "PROS" = 80.5, 
                     "SEMA3" = 83.5, 
                     "SELE" = 86, 
                     "SEMA6" = 87.5, 
                     "SEMA7" = 89, 
                     "TWEAK" = 90, 
                     "VEGF" = 92.5) 
gene_grodown <- ggplot(gc.s.rdown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown

p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.1 EC的source下调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rdown.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5, 25.5, 27.5, 28.5,
                          29.5, 31.5, 33.5, 34.5, 37.5, 41.5, 42.5, 67.5, 68.5,
                          69.5, 71.5, 72.5, 77.5, 79.5, 81.5, 85.5, 86.5, 88.5,
                          89.5, 90.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC的source下调受体气泡图B


####EC的下调受体气泡图（target）#####################
#EC的下调受体气泡图（target）
gc.t.rdown <- receptor.down %>%  
  filter(target == "EC")

gc.t.rdown <- gc.t.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.rdown$pathway_name)
gc.t.rdown.use = gc.t.rdown[, "interaction_name", drop = F]
gc.t.rdown$interaction_name <- as.factor(gc.t.rdown$interaction_name)
gc.t.rdown$interaction_name <- fct_inorder(gc.t.rdown$interaction_name)
gc.t.rdown$interaction_name_2 <- as.factor(gc.t.rdown$interaction_name_2)
gc.t.rdown$interaction_name_2 <- fct_inorder(gc.t.rdown$interaction_name_2)
gc.t.rdown$pathway_name <- as.factor(gc.t.rdown$pathway_name)
gc.t.rdown$pathway_name <- fct_inorder(gc.t.rdown$pathway_name)
gc.t.rdown <- gc.t.rdown %>% 
  mutate(p="")

gc.t.rdown$signway <- paste(gc.t.rdown$source,
                          "  ->  ", 
                          gc.t.rdown$target) 

table(unique(gc.t.rdown$interaction_name_2))
summary_data <- gc.t.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.t.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Receptor of EC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(1.5, 2.5, 5.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, 
                     "SEMA3" = 2, 
                     "TGFb" = 4, 
                     "VEGF" = 10)  
gene_grodown <- ggplot(gc.t.rdown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown
p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.3 EC-GC的target下调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rdown.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 5.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-GC的target下调受体气泡图B


