rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(clustree))install.packages("clustree")
  if(!require(devtools))install.packages("devtools")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(Matrix))install.packages("Matrix")
  if(!require(SeuratDisk))remotes::install_github("mojaveazure/seurat-disk")
  if(!require(Rcpp))install.packages("Rcpp")
  if(!require(RcppEigen))install.packages("RcppEigen")
  if(!require(ggsci))install.packages("ggsci")
  if(!require(viridis))install.packages("viridis")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(scibetR))devtools::install_github("zwj-tina/scibetR")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(MySeuratWrappers))install.packages("MySeuratWrappers")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/确定注释.rds")
Idents(mouse_all) <- mouse_all$usetype



####设定marker######################
#设定marker
EC.marker <- c("Egf17", "Cdh5", "Ptprb", "Flt1", "Pecam1", "Cd300lg", 
               "Cyyr1", "Rasgrp3", "F8", "Gpihbp1", "Kdr", "Emcn", 
               "Ehd3", "Srgn", "Egfl7", "S100a10", "Clic4", "Ly6a")
IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
PTC.marker <- c("Tacstd2", "Pcbd1", "Cdh16", "Pdzk1ip1", "Atp1b1", "Spp1", 
                "Fxyd2", "Cryab", "Ldhb", "S100g", "Slc12a3")
POD.marker <- c("Nphs1", "Enpep", "Thsd7a", "Dpp4", "Npnt", "Plce1", 
                "Ptpro", "Nphs2", "Clic3", "Cdkn1c", "Synpo", "Mafb", 
                "Wt1", "Cd2ap", "Podxl", "Col4a3", "Golim4", "Magi2", 
                "Pard3b")
MC.marker <- c("Gucy1a3", "Gata3", "Pdgfrb", "Sfrp2", "Itga8", "Myl9", 
               "Agtr1a", "Hopx", "Ptn", "Mgp", "Ren1", "Tpm2",
               "Cald1", "Nr4a2", "Actn1", "Akap12", "Rock1", "Sh3bgrl")
FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")

ALL.marker <- unique(c(EC.marker, IMMUNE.marker, PTC.marker, POD.marker, 
                       MC.marker, FIB.marker, VSMC.marker))
  

####UMAP图###############################
#UMAP图
cell_colors <- c(  
  "Immune Cell" = "#E6194B", # 红色  
  "EC" = "#3CB44B",         # 绿色  
  "PTC" = "#FFE119",        # 黄色  
  "POD" = "#4363D8",        # 蓝色  
  "MC" = "#F58231",         # 橙色  
  "FIB" = "#911EBB",        # 紫色  
  "VSMC" = "navy"        # 青色  
)  
DimPlot(mouse_all, reduction = "umap", group.by = "usetype", 
        cols = cell_colors,   
        label = TRUE, pt.size = 0.5) +
  NoLegend()
  ##得到1.1 UMAP图-确定注释 

DimPlot(mouse_all, reduction = "umap", group.by = "usetype", 
        cols = cell_colors,   
        label = FALSE, pt.size = 0.5)
  ##得到1.2 UMAP图-确定注释（图例）

FeaturePlot(mouse_all, features = EC.marker, raster=FALSE, ncol = 4)
  ##得到2.UMAP图-EC.marker

FeaturePlot(mouse_all, features = IMMUNE.marker, raster=FALSE, ncol = 4)
  ##得到3.UMAP图-IMMUNE.marker

FeaturePlot(mouse_all, features = PTC.marker, raster=FALSE, ncol = 4)
  ##得到4.UMAP图-PTC.marker

FeaturePlot(mouse_all, features = POD.marker, raster=FALSE, ncol = 4)
  ##得到5.UMAP图-POD.marker

FeaturePlot(mouse_all, features = MC.marker, raster=FALSE, ncol = 4)
  ##得到6.UMAP图-MC.marker

FeaturePlot(mouse_all, features = FIB.marker, raster=FALSE, ncol = 4)
  ##得到7.UMAP图-FIB.marker

FeaturePlot(mouse_all, features = VSMC.marker, raster=FALSE, ncol = 4)
  ##得到8.UMAP图-VSMC.marker


####点图和小提琴图##################################
#点图和小提琴图
DotPlot(mouse_all, features = ALL.marker) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = c(18.5, 28.5, 39.5, 57.5, 74.5, 78.5))  
  ##得到9.点图-全部marker

VlnPlot(mouse_all, 
        features = ALL.marker,
        stacked=T,
        pt.size=0,
        direction = "horizontal", 
        x.lab = '', y.lab = '')+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
  ##得到10.小提琴图-全部marker



####spearman相关性##########################
#spearman相关性
table(mouse_all$usetype)
mouse_pearson <- AverageExpression(mouse_all,
                                   group.by = "usetype",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),1000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到11.全部细胞注释的spearman相关性
