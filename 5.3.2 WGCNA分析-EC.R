rm(list = ls())
gc()

#载入必要的package
{
  if(!require(BiocManager))install.packages("BiocManager")
  if(!require(WGCNA))BiocManager::install("WGCNA")
  if(!require(igraph))BiocManager::install("igraph")
  if(!require(devtools))BiocManager::install("devtools")
  if(!require(GeneOverlap))BiocManager::install("GeneOverlap")
  if(!require(ggrepel))BiocManager::install("ggrepel")
  if(!require(UCell))BiocManager::install("UCell")
  if(!require(ggforestplot))devtools::install_github("NightingaleHealth/ggforestplot")
  if(!require(hdWGCNA))devtools::install_github('smorabit/hdWGCNA', ref='dev')
  if(!require(ggrastr))install.packages("ggrastr")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(cowplot))install.packages("cowplot")
  if(!require(patchwork))install.packages("patchwork")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(stringr))install.packages("stringr")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")
mouse_all <- subset(mouse_all, usesubtype == "EC")
mouse_all$usesubtype <- droplevels(mouse_all$usesubtype,
                                   exclude = setdiff(
                                     levels(mouse_all$usesubtype),
                                     unique(mouse_all$usesubtype)))

table(mouse_all$usesubtype)
Idents(mouse_all) <- mouse_all$usesubtype



####环境准备#########################
#环境准备
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)


DefaultAssay(mouse_all) <- 'RNA'
DimPlot(mouse_all, label = T)
DimPlot(mouse_all, group.by = "orig.ident", label = T)



####WGCNA对象的预处理##########################
#WGCNA对象的预处理
mouse_all <- SetupForWGCNA(mouse_all,
                           wgcna_name = "Mouse",
                           gene_select = "fraction",
                           fraction = 0.05)

dim(mouse_all)

mouse_all <- MetacellsByGroups(seurat_obj =  mouse_all,
                               group.by = c("usesubtype", "group"),
                               k = 10,
                               max_shared = 20,
                               reduction = "umap",
                               ident.group = "usesubtype",
                               min_cells = 10)

mouse_all <- NormalizeMetacells(mouse_all)



####选择软阈值################################
#选择软阈值
mouse_all <- SetDatExpr(mouse_all, assay = 'RNA', slot = 'data')

mouse_all <- TestSoftPowers(mouse_all, networkType = 'signed')
plot <- PlotSoftPowers(mouse_all)
wrap_plots(plot, ncol=2)
  ##得到1.EC亚群软阈值折线图

mouse_all <- ConstructNetwork(mouse_all,
                              soft_power = 12,
                              tom_name = "Mouse_Test",
                              setDatExpr = F)

PlotDendrogram(mouse_all, main='scRNA hdWGCNA Dendrogram')
  ##得到2.EC亚群模块树状图

TOM <- GetTOM(mouse_all)



####计算特征基因和连通性#######################
#计算特征基因和连通性
mouse_all <- ScaleData(mouse_all)
mouse_all <- ModuleEigengenes(mouse_all, group.by.vars = 'group')
hMEs <- GetMEs(mouse_all)
MEs <- GetMEs(mouse_all, harmonized=FALSE)
mouse_all <- ModuleConnectivity(mouse_all)
modules <- GetModules(mouse_all)
write.csv(modules, "EC亚群基因模块列表.csv")
  ##得到EC亚群基因模块列表


##重命名module
mouse_all <- ResetModuleNames(mouse_all,new_name = "M")


##修改module颜色
mod_color_df <- GetModules(mouse_all) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)
  ##得到仅包含模块及其唯一颜色的表格

n_mods <- nrow(mod_color_df) - 1
  ##删除第一个模块

newcolor <- c("black", "blue", "greenyellow", "darkorange", "grey60", 
              "skyblue", "tan", "brown", "royalblue", "red", "purple", 
              "yellow", "lightyellow", "cyan", "turquoise", "orange", 
              "pink", "green", "lightgreen", "darkturquoise", "darkgreen",
              "lightcyan", "salmon", "midnightblue", "darkred", "magenta",
              "darkgrey", "#FF6347")
mouse_all <- ResetModuleColors(mouse_all, newcolor)

modules <- GetModules(mouse_all)
  ##得到新的模块

hub_25 <- GetHubGenes(mouse_all, n_hubs = 25)
write.csv(hub_25, "EC亚群hub基因-top25.csv")
  ##得到EC亚群hub基因-top25

hub_50 <- GetHubGenes(mouse_all, n_hubs = 50)
write.csv(hub_50, "EC亚群hub基因-top50.csv")
  ##得到EC亚群hub基因-top50

hub_all <- GetHubGenes(mouse_all, n_hubs = 1000)
write.csv(hub_all, "EC亚群hub基因-all.csv")
  ##得到EC亚群hub基因-all


##使用UCell计算TOP50的hub基因
mouse_all <- ModuleExprScore(mouse_all,
                             n_genes = 50,
                             method = "UCell")



####网络可视化################################
#网络可视化
PlotKMEs(mouse_all, ncol = 8)
  ##得到3.EC亚群模块基因的kME排序图


plot_hMEs <- ModuleFeaturePlot(mouse_all,
                               reduction = "umap",
                               features = "hMEs", 
                               order = TRUE,
                               raster = T)
wrap_plots(plot_hMEs, ncol = 8)
  ##得到4.EC亚群协调模块的特征基因的UMAP图


plot_score <- ModuleFeaturePlot(mouse_all,
                                reduction = "umap",
                                features = "hMEs",
                                order = TRUE,
                                raster = T,
                                ucell = TRUE)
wrap_plots(plot_score, ncol = 8)
  ##得到5.EC亚群Ucell的特征基因的UMAP图


mouse_all@meta.data <- cbind(mouse_all@meta.data,
                             GetMEs(mouse_all, harmonized=TRUE))
MEs <- GetMEs(mouse_all, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- as.numeric(gsub("M(\\d+)", "\\1", mods))  
mods <- mods[order(mods_num)]  

DotPlot(mouse_all,
        features = mods,
        group.by = 'metacell_grouping') + 
  coord_flip()+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5))
  ##得到7.模块在EC亚群表达的的点图


p1 <- VlnPlot(mouse_all,
              features = mods,
              group.by = 'metacell_grouping',
              pt.size = 0,
              ncol = 8)
p1= p1+geom_boxplot(width=.25, fill='white')
p1 <- p1 + xlab('') + ylab('hME') + NoLegend()
p1
  ##得到8.模块在EC亚群表达的的小提琴图


ModuleNetworkPlot(mouse_all,
                  outdir = "./EC亚群各模块网状图")
  ##在文件夹得到EC亚群各模块的网状图
