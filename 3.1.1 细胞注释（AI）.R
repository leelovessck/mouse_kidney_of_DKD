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
  if(!require(GPTCelltype))remotes::install_github("Winnie09/GPTCelltype")
  if(!require(openai))install.packages("openai")
}

setwd("E:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/样本分群.rds")
Idents(mouse_all) <- mouse_all$seurat_clusters



####细胞marker####################
#细胞marker
markers.all <- FindAllMarkers(mouse_all,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)



####GPT调用#################
#GPT调用
Sys.setenv(OPENAI_API_KEY = '')
res <- gptcelltype(markers.all, 
                   tissuename = 'mouse glomerulus', 
                   model = 'gpt-4',
                   topgenenumber = 10)
res
##"Identify cell types of human PBMC cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. \n0:CD7,GNLY,CCL5\n1:S100A8,TYMP,S100A9\n2:HLA-DPB1,MS4A1,HLA-DQB1"



####细胞注释（AI,TOP20）##########################
#AIcelltype(TOP20)

##cluster0代表壁层上皮细胞或其他间质细胞：
##Sapcd1等基因可能与这些细胞的特性或功能相关。

##cluster1代表足细胞：
##Itga7和Gucy1a2等基因的表达在足细胞中较为常见。

##cluster2代表系膜细胞

##cluster3肾小球内皮细胞：
##Clec4e等基因与内皮细胞的功能和特性紧密相关


##cluster4代表上皮细胞或其他类型的肾小球细胞

##cluster5代表足细胞：
##Cdh1（E-cadherin）是上皮细胞的一个典型标记，特别是在足细胞中表达较高。

##cluster6代表上皮细胞或其他类型的肾小球细胞

##cluster7代表与免疫相关的细胞，如巨噬细胞或树突状细胞：
##C3ar1、Siglec1等基因与免疫细胞的功能相关。

##cluster8代表其他类型的肾小球辅助细胞或混合细胞类型：

##cluster9代表免疫细胞，特别是B淋巴细胞或浆细胞：
##Iglv1、Hba-a2等基因明显与免疫细胞相关

##cluster10代表发育中的肾小球细胞或其他非特定类型的细胞

AIcelltype <- c("0" = "PTC",
                "1" = "POD",
                "2" = "MC",
                "3" = "GEC",
                "4" = "Mixed",
                "5" = "POD",
                "6" = "Mixed",
                "7" = "Immune Cell",
                "8" = "Mixed",
                "9" = "Immune Cell",
                "10" = "Growing Cell")
names(AIcelltype) <- levels(mouse_all)
mouse_all <- RenameIdents(mouse_all, AIcelltype)
mouse_all@meta.data$AIcelltype <- Idents(mouse_all)
Idents(mouse_all) <- mouse_all$AIcelltype



####可视化（AI,TOP20）##########################
#可视化（AI,TOP20）
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "green", "cyan",
                 "yellow", "purple", "orange"))
  ##得到1.UMAP-AI注释（TOP20）

table(mouse_all$AIcelltype)
mouse_pearson <- AverageExpression(mouse_all,
                                   group.by = "AIcelltype",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),20000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到2.AI注释的spearman相关性（TOP20）

EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mouse_all, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.AI注释内皮细胞marker（TOP20）

MC.marker <- c("Gata3", "Pdgfrb", "Myl9", "Acta2", "Ren1")
p1 <- FeaturePlot(mouse_all, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到4.AI注释系膜细胞marker（TOP20）

POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
p1 <- FeaturePlot(mouse_all, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = POD.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到5.AI注释足细胞marker（TOP20）

IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mouse_all, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到6.AI注释免疫细胞marker（TOP20）

PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
p1 <- FeaturePlot(mouse_all, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到7.AI注释壁上皮细胞marker（TOP20）

VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")
p1 <- FeaturePlot(mouse_all, features = VSMC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = VSMC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到8.AI注释血管平滑肌细胞marker（TOP20）

FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
p1 <- FeaturePlot(mouse_all, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = FIB.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到9.AI注释成纤维细胞marker（TOP20）
