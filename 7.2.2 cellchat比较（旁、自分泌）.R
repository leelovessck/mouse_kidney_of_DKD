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
}



####数据准备############################
#数据准备
con.chat <- readRDS("./data source/cellchat/control（自、旁分泌）.rds")
dkd.chat <- readRDS("./data source/cellchat/DKD（自、旁分泌）.rds")
con.chat <- netAnalysis_computeCentrality(con.chat) 
dkd.chat <- netAnalysis_computeCentrality(dkd.chat) 
object.list <- list(CON = con.chat, DKD = dkd.chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



####柱状图#############################
#柱状图
##基础柱状图
gg1 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2))
gg2 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2), 
                           measure = "weight")
gg1 + gg2
  ##得到1.旁、自分泌比较柱状图（结论）



##对比柱状图
gg1 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "count",
               stacked = T, 
               do.stat = TRUE)
gg2 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "count",
               stacked = F, 
               do.stat = TRUE)
gg1 + gg2
  ##得到2.1 旁、自分泌的比较柱状图-数量（详细通路）

gg1 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "count",
               stacked = F, 
               do.stat = TRUE,
               sources.use = "EC")
gg2 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "count",
               stacked = F, 
               do.stat = TRUE,
               targets.use = "EC")
gg1 + gg2
  ##得到2.2 EC的旁、自分泌的比较柱状图-数量（详细通路）


gg1 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "weight",
               stacked = T, 
               do.stat = TRUE)
gg2 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "weight",
               stacked = F, 
               do.stat = TRUE)
gg1 + gg2
  ##得到3.1 旁、自分泌的比较柱状图-强度（详细通路）

gg1 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "weight",
               stacked = F, 
               do.stat = TRUE,
               sources.use = "EC")
gg2 <- rankNet(cellchat,
               mode = "comparison", 
               measure = "weight",
               stacked = F, 
               do.stat = TRUE,
               targets.use = "EC")
gg1 + gg2
  ##得到3.2 EC的旁、自分泌的比较柱状图-强度（详细通路）



####环形图#############################
#环形图
##直接显示两组相减
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, 
                          weight.scale = T)
netVisual_diffInteraction(cellchat, 
                          weight.scale = T, 
                          measure = "weight")
  ##得到4.1 组间旁、自分泌环形图

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, 
                          weight.scale = T,
                          sources.use = "EC")
netVisual_diffInteraction(cellchat, 
                          weight.scale = T, 
                          measure = "weight",
                          targets.use = "EC")
  ##得到4.2 组间旁、自分泌环形图(EC)



####热图#########################
#热图
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
  ##得到5.旁、自分泌的比较热图



####细胞间通讯强度####################
#细胞间通讯强度
levels(object.list[[1]]@idents) 
levels(object.list[[2]]@idents) 
group.cellType <- levels(object.list[[1]]@idents)
group.cellType <- factor(group.cellType, levels = levels(object.list[[1]]@idents))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})#合并细胞类型


##展示通讯强度
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
  ##得到6.旁、自分泌的in out散点图


##展示通讯差异
for (i in seq_along(group.cellType)) {  
  ident <- group.cellType[i]  
  assign(paste0("gg", i), netAnalysis_signalingChanges_scatter(cellchat, idents.use = ident), envir = .GlobalEnv)  
}
p1 <- patchwork::wrap_plots(plots = list(gg1, gg2, gg3, gg4, gg5, gg6, gg7),
                            ncol = 4)
filename <- "7.全部通路的in out散点图(按细胞类型).png"  
png(filename, width = 2000, height = 750)  
print(p1)  
dev.off()



##outgoing热图
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, 
                       object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 10, 
                                        height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 10, 
                                        height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  ##得到8.1 旁、自分泌的out热图(按细胞类型)


##incoming热图
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 10, 
                                        height = 15, 
                                        color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "incoming",
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1],
                                        width = 10,
                                        height = 15, 
                                        color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  ##得到8.2 旁、自分泌的in热图(按细胞类型)
