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
write.csv(net, "改变受配体.csv")



####TGFb通路###############################
#TGFb通路
net_TGFb <- net %>%  
  filter(pathway_name == "TGFb" & (source %in% "EC" | target %in% "EC")) 

unique(net_TGFb$pathway_name)
net_TGFb.use = net_TGFb[, "interaction_name", drop = F]

net_TGFb$interaction_name <- as.factor(net_TGFb$interaction_name)
net_TGFb$interaction_name <- fct_inorder(net_TGFb$interaction_name)
net_TGFb$interaction_name_2 <- as.factor(net_TGFb$interaction_name_2)
net_TGFb$interaction_name_2 <- fct_inorder(net_TGFb$interaction_name_2)
net_TGFb$pathway_name <- as.factor(net_TGFb$pathway_name)
net_TGFb$pathway_name <- fct_inorder(net_TGFb$pathway_name)
net_TGFb <- net_TGFb %>% 
  mutate(p="")
net_TGFb$signway <- paste(net_TGFb$source,
                          "  ->  ", 
                          net_TGFb$target) 

p1 <- ggplot(net_TGFb,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "TGFb pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到1.1 EC的TGFb气泡图A

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_TGFb.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "TGFb pathway")
gg2
  ##得到1.2 EC（target）的TGFb气泡图B



####PTN通路###############################
#PTN通路
net_PTN <- net %>%  
  filter(pathway_name == "PTN" & (source %in% "EC" | target %in% "EC")) 

unique(net_PTN$pathway_name)
net_PTN.use = net_PTN[, "interaction_name", drop = F]

net_PTN$interaction_name <- as.factor(net_PTN$interaction_name)
net_PTN$interaction_name <- fct_inorder(net_PTN$interaction_name)
net_PTN$interaction_name_2 <- as.factor(net_PTN$interaction_name_2)
net_PTN$interaction_name_2 <- fct_inorder(net_PTN$interaction_name_2)
net_PTN$pathway_name <- as.factor(net_PTN$pathway_name)
net_PTN$pathway_name <- fct_inorder(net_PTN$pathway_name)
net_PTN <- net_PTN %>% 
  mutate(p="")
net_PTN$signway <- paste(net_PTN$source,
                         "  ->  ", 
                         net_PTN$target) 

p1 <- ggplot(net_PTN,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "PTN pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到2.1 EC的PTN气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_PTN.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "PTN pathway")
gg1
  ##得到2.2 EC（source）的PTN气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_PTN.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "PTN pathway")
gg2
  ##得到2.3 EC（target）的PTN气泡图B



####CSF通路###############################
#CSF通路
net_CSF <- net %>%  
  filter(pathway_name == "CSF" & (source %in% "EC" | target %in% "EC")) 

unique(net_CSF$pathway_name)
net_CSF.use = net_CSF[, "interaction_name", drop = F]

net_CSF$interaction_name <- as.factor(net_CSF$interaction_name)
net_CSF$interaction_name <- fct_inorder(net_CSF$interaction_name)
net_CSF$interaction_name_2 <- as.factor(net_CSF$interaction_name_2)
net_CSF$interaction_name_2 <- fct_inorder(net_CSF$interaction_name_2)
net_CSF$pathway_name <- as.factor(net_CSF$pathway_name)
net_CSF$pathway_name <- fct_inorder(net_CSF$pathway_name)
net_CSF <- net_CSF %>% 
  mutate(p="")
net_CSF$signway <- paste(net_CSF$source,
                         "  ->  ", 
                         net_CSF$target) 

p1 <- ggplot(net_CSF,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "CSF pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到3.1 EC的CSF气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_CSF.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "CSF pathway")
gg1
  ##得到3.2 EC（source）的CSF气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_CSF.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "CSF pathway")
gg2
  ##得到3.3 EC（target）的CSF气泡图B



####EGF通路###############################
#EGF通路
net_EGF <- net %>%  
  filter(pathway_name == "EGF" & (source %in% "EC" | target %in% "EC")) 

unique(net_EGF$pathway_name)
net_EGF.use = net_EGF[, "interaction_name", drop = F]

net_EGF$interaction_name <- as.factor(net_EGF$interaction_name)
net_EGF$interaction_name <- fct_inorder(net_EGF$interaction_name)
net_EGF$interaction_name_2 <- as.factor(net_EGF$interaction_name_2)
net_EGF$interaction_name_2 <- fct_inorder(net_EGF$interaction_name_2)
net_EGF$pathway_name <- as.factor(net_EGF$pathway_name)
net_EGF$pathway_name <- fct_inorder(net_EGF$pathway_name)
net_EGF <- net_EGF %>% 
  mutate(p="")
net_EGF$signway <- paste(net_EGF$source,
                         "  ->  ", 
                         net_EGF$target) 

p1 <- ggplot(net_EGF,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "EGF pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到4.1 EC的EGF气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_EGF.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "EGF pathway")
gg1
  ##得到4.2 EC（source）的EGF气泡图B-无

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_EGF.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "EGF pathway")
gg2
  ##得到4.3 EC（target）的EGF气泡图B



####FGF通路###############################
#FGF通路
net_FGF <- net %>%  
  filter(pathway_name == "FGF" & (source %in% "EC" | target %in% "EC")) 

unique(net_FGF$pathway_name)
net_FGF.use = net_FGF[, "interaction_name", drop = F]

net_FGF$interaction_name <- as.factor(net_FGF$interaction_name)
net_FGF$interaction_name <- fct_inorder(net_FGF$interaction_name)
net_FGF$interaction_name_2 <- as.factor(net_FGF$interaction_name_2)
net_FGF$interaction_name_2 <- fct_inorder(net_FGF$interaction_name_2)
net_FGF$pathway_name <- as.factor(net_FGF$pathway_name)
net_FGF$pathway_name <- fct_inorder(net_FGF$pathway_name)
net_FGF <- net_FGF %>% 
  mutate(p="")
net_FGF$signway <- paste(net_FGF$source,
                         "  ->  ", 
                         net_FGF$target) 

p1 <- ggplot(net_FGF,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "FGF pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到5.1 EC的FGF气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_FGF.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "FGF pathway")
gg1
  ##得到5.2 EC（source）的FGF气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_FGF.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "FGF pathway")
gg2
  ##得到5.3 EC（target）的FGF气泡图B



####IGF通路###############################
#IGF通路
net_IGF <- net %>%  
  filter(pathway_name == "IGF" & (source %in% "EC" | target %in% "EC")) 

unique(net_IGF$pathway_name)
net_IGF.use = net_IGF[, "interaction_name", drop = F]

net_IGF$interaction_name <- as.factor(net_IGF$interaction_name)
net_IGF$interaction_name <- fct_inorder(net_IGF$interaction_name)
net_IGF$interaction_name_2 <- as.factor(net_IGF$interaction_name_2)
net_IGF$interaction_name_2 <- fct_inorder(net_IGF$interaction_name_2)
net_IGF$pathway_name <- as.factor(net_IGF$pathway_name)
net_IGF$pathway_name <- fct_inorder(net_IGF$pathway_name)
net_IGF <- net_IGF %>% 
  mutate(p="")
net_IGF$signway <- paste(net_IGF$source,
                         "  ->  ", 
                         net_IGF$target) 

p1 <- ggplot(net_IGF,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "IGF pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到6.1 EC的IGF气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_IGF.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "IGF pathway")
gg1
  ##得到6.2 EC（source）的IGF气泡图B-无

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_IGF.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "IGF pathway")
gg2
  ##得到6.3 EC（target）的IGF气泡图B



####VEGF通路###############################
#VEGF通路
net_VEGF <- net %>%  
  filter(pathway_name == "VEGF" & (source %in% "EC" | target %in% "EC")) 

unique(net_VEGF$pathway_name)
net_VEGF.use = net_VEGF[, "interaction_name", drop = F]

net_VEGF$interaction_name <- as.factor(net_VEGF$interaction_name)
net_VEGF$interaction_name <- fct_inorder(net_VEGF$interaction_name)
net_VEGF$interaction_name_2 <- as.factor(net_VEGF$interaction_name_2)
net_VEGF$interaction_name_2 <- fct_inorder(net_VEGF$interaction_name_2)
net_VEGF$pathway_name <- as.factor(net_VEGF$pathway_name)
net_VEGF$pathway_name <- fct_inorder(net_VEGF$pathway_name)
net_VEGF <- net_VEGF %>% 
  mutate(p="")
net_VEGF$signway <- paste(net_VEGF$source,
                         "  ->  ", 
                         net_VEGF$target) 

p1 <- ggplot(net_VEGF,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "VEGF pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到7.1 EC的VEGF气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_VEGF.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "VEGF pathway")
gg1
  ##得到7.2 EC（source）的VEGF气泡图B-无

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_VEGF.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "VEGF pathway")
gg2
  ##得到7.3 EC（target）的VEGF气泡图B



####COLLAGEN通路###############################
#COLLAGEN通路
net_COLLAGEN <- net %>%  
  filter(pathway_name == "COLLAGEN" & (source %in% "EC" | target %in% "EC")) 

unique(net_COLLAGEN$pathway_name)
net_COLLAGEN.use = net_COLLAGEN[, "interaction_name", drop = F]

net_COLLAGEN$interaction_name <- as.factor(net_COLLAGEN$interaction_name)
net_COLLAGEN$interaction_name <- fct_inorder(net_COLLAGEN$interaction_name)
net_COLLAGEN$interaction_name_2 <- as.factor(net_COLLAGEN$interaction_name_2)
net_COLLAGEN$interaction_name_2 <- fct_inorder(net_COLLAGEN$interaction_name_2)
net_COLLAGEN$pathway_name <- as.factor(net_COLLAGEN$pathway_name)
net_COLLAGEN$pathway_name <- fct_inorder(net_COLLAGEN$pathway_name)
net_COLLAGEN <- net_COLLAGEN %>% 
  mutate(p="")
net_COLLAGEN$signway <- paste(net_COLLAGEN$source,
                          "  ->  ", 
                          net_COLLAGEN$target) 

p1 <- ggplot(net_COLLAGEN,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "COLLAGEN pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到8.1 EC的COLLAGEN气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_COLLAGEN.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "COLLAGEN pathway")
gg1
  ##得到8.2 EC（source）的COLLAGEN气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_COLLAGEN.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "COLLAGEN pathway")
gg2
  ##得到8.3 EC（target）的COLLAGEN气泡图B



####LAMININ通路###############################
#LAMININ通路
net_LAMININ <- net %>%  
  filter(pathway_name == "LAMININ" & (source %in% "EC" | target %in% "EC")) 

unique(net_LAMININ$pathway_name)
net_LAMININ.use = net_LAMININ[, "interaction_name", drop = F]

net_LAMININ$interaction_name <- as.factor(net_LAMININ$interaction_name)
net_LAMININ$interaction_name <- fct_inorder(net_LAMININ$interaction_name)
net_LAMININ$interaction_name_2 <- as.factor(net_LAMININ$interaction_name_2)
net_LAMININ$interaction_name_2 <- fct_inorder(net_LAMININ$interaction_name_2)
net_LAMININ$pathway_name <- as.factor(net_LAMININ$pathway_name)
net_LAMININ$pathway_name <- fct_inorder(net_LAMININ$pathway_name)
net_LAMININ <- net_LAMININ %>% 
  mutate(p="")
net_LAMININ$signway <- paste(net_LAMININ$source,
                              "  ->  ", 
                              net_LAMININ$target) 

p1 <- ggplot(net_LAMININ,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "LAMININ pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到9.1 EC的LAMININ气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_LAMININ.use, 
                        sources.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "LAMININ pathway")
gg1
  ##得到9.2 EC（source）的LAMININ气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_LAMININ.use, 
                        targets.use = "EC",
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "LAMININ pathway")
gg2
  ##得到9.3 EC（target）的LAMININ气泡图B
