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

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/确定注释.rds")
Idents(mouse_all) <- mouse_all$usetype



####免疫细胞降维分群########################
#免疫细胞降维分群
sub_immune <- subset(mouse_all, usetype == "Immune Cell")
sub_immune <- FindClusters(sub_immune, resolution = 0.5)
sub_immune <- RunUMAP(sub_immune, dims = 1:25)



####marker可视化##########################
#marker可视化
Macrophage.marker <- c("Stab1", "Apoe", "Lyz2", "Tgfbi", "Vcan")
p1 <- FeaturePlot(sub_immune, features = Macrophage.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Macrophage.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.巨噬细胞marker

B.marker <- c("Pax5", "Ebf1", "Cd79a", "Cd79b")
p1 <- FeaturePlot(sub_immune, features = B.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = B.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = B.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.B细胞marker

T.marker <- c("Skap1", "Cd247", "Cd3g", "Cd8b1", "Cd8a", "Cd4")
p1 <- FeaturePlot(sub_immune, features = T.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = T.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.T细胞marker

DC.marker <- c("Clec9a", "Ppt1", "Itgam", "Sirpa", "Itgax", "H2-Aa", 
               "H2-Ab1", "Cd74", "Ccl5")
p1 <- FeaturePlot(sub_immune, features = DC.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = DC.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = DC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到4.DC细胞marker

Erythrocyte.marker <- c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")
p1 <- FeaturePlot(sub_immune, features = Erythrocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Erythrocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Erythrocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到5.AI注释红细胞marker

Granulocyte.marker <- c("Thbs1", "Cd14", "Junb", "Lcn2", "S100a9", 
                        "Il1b", "Ly6g")
p1 <- FeaturePlot(sub_immune, features = Granulocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Granulocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Granulocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到6.粒细胞marker

NK.marker <- c("Nkg7", "Gzmb", "Gzma")
p1 <- FeaturePlot(sub_immune, features = NK.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = NK.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = NK.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到7.NK细胞marker



####细胞注释（hand）##########################
#handcelltype
handcelltype <- c("0" = "Mixed（MAC/T）",
                  "1" = "Mixed（DC/T/G）",
                  "2" = "Erythrocyte",
                  "3" = "Mixed（Erythrocyte/B）",
                  "4" = "Mixed（MAC/G）")
names(handcelltype) <- levels(sub_immune)
sub_immune <- RenameIdents(sub_immune, handcelltype)
sub_immune@meta.data$handcelltype <- Idents(sub_immune)
Idents(sub_immune) <- sub_immune$handcelltype
saveRDS(sub_immune, "./data source/rds/免疫细胞（含mixed）.rds")



####可视化（hand）##########################
#可视化（hand）
DimPlot(sub_immune,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "green", "cyan", "purple"))
  ##得到8.UMAP-手动注释

table(sub_immune$handcelltype)
mouse_pearson <- AverageExpression(sub_immune,
                                   group.by = "handcelltype",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),20000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到9.手动注释的spearman相关性


p1 <- FeaturePlot(sub_immune, features = Macrophage.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Macrophage.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.1 手动注释巨噬细胞marker

p1 <- FeaturePlot(sub_immune, features = B.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = B.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = B.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.2 手动注释B细胞marker

p1 <- FeaturePlot(sub_immune, features = T.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = T.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.3 手动注释T细胞marker

p1 <- FeaturePlot(sub_immune, features = DC.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = DC.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = DC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.4 手动注释DC细胞marker

p1 <- FeaturePlot(sub_immune, features = Erythrocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Erythrocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Erythrocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.5 手动注释红细胞marker

Granulocyte.marker <- c("Thbs1", "Cd14", "Junb", "Lcn2", "S100a9", 
                        "Il1b", "Ly6g")
p1 <- FeaturePlot(sub_immune, features = Granulocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Granulocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Granulocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.6 手动注释粒细胞marker
