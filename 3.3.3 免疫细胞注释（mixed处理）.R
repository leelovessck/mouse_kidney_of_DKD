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
mouse_all <- readRDS("./data source/rds/免疫细胞（含mixed）.rds")
Idents(mouse_all) <- mouse_all$handcelltype

####设定marker##############################
#设定marker
Macrophage.marker <- c("Stab1", "Apoe", "Lyz2", "Tgfbi", "Vcan")
B.marker <- c("Pax5", "Ebf1", "Cd79a", "Cd79b")
T.marker <- c("Skap1", "Cd247", "Cd3g", "Cd8b1", "Cd8a", "Cd4")
DC.marker <- c("Clec9a", "Ppt1", "Itgam", "Sirpa", "Itgax", "H2-Aa", 
               "H2-Ab1", "Cd74", "Ccl5")
Erythrocyte.marker <- c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")
Granulocyte.marker <- c("Thbs1", "Cd14", "Junb", "Lcn2", "S100a9", 
                        "Il1b", "Ly6g")



####筛选cluster0##########################
#筛选cluster0
mixed_1 <- subset(mouse_all, handcelltype == "Mixed（MAC/T）")
object_for_clustree <- mixed_1
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到1.1 Mixed（cluster0）分群树状图

mixed_1 <- FindClusters(mixed_1, resolution = 1.2)
mixed_1 <- RunUMAP(mixed_1, dims = 1:25)

DimPlot(mixed_1,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到1.2 UMAP-Mixed（cluster0）分群



####cluster0的marker可视化##########################
#cluster0的marker可视化
p1 <- FeaturePlot(mixed_1, features = Macrophage.marker)
p2 <- DotPlot(mixed_1, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = Macrophage.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到1.3 巨噬细胞细胞marker

p1 <- FeaturePlot(mixed_1, features = T.marker)
p2 <- DotPlot(mixed_1, features = T.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到1.4 T细胞marker



####细胞注释Mixed（cluster0）##########################
#handcelltype
subtype <- c("0" = "Macrophage",
             "1" = "EC",
             "2" = "EC",
             "3" = "Macrophage",
             "4" = "Macrophage",
             "5" = "Macrophage")
names(subtype) <- levels(mixed_1)
mixed_1 <- RenameIdents(mixed_1, subtype)
mixed_1@meta.data$subtype <- Idents(mixed_1)
Idents(mixed_1) <- mixed_1$subtype


p1 <- FeaturePlot(mixed_1, features = c(Macrophage.marker, T.marker),
                  raster=FALSE)
p2 <- DotPlot(mixed_1, features = c(Macrophage.marker, T.marker))+
  coord_flip()
p3 <- VlnPlot(mixed_1, features = c(Macrophage.marker, T.marker),
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到1.5 Mixed（cluster0）的marker

Idents(mouse_all, cells = colnames(mixed_1)) <- Idents(mixed_1)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()



####筛选Mixed（cluster1）##########################
#筛选Mixed（cluster1）
mixed_2 <- subset(mouse_all, handcelltype == "Mixed（DC/T/G）")
object_for_clustree <- mixed_2
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到2.1 Mixed（cluster1）分群树状图

mixed_2 <- FindClusters(mixed_2, resolution = 1.2)
mixed_2 <- RunUMAP(mixed_2, dims = 1:25)

DimPlot(mixed_2,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到2.2 UMAP-Mixed（cluster1）分群



####cluster2的marker可视化##########################
#cluster2的marker可视化
p1 <- FeaturePlot(mixed_2, features = Macrophage.marker)
p2 <- DotPlot(mixed_2, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = Macrophage.marker,
              pt.size = 0, log = TRUE )
p1 / p2 / p3
  ##得到2.3 巨噬细胞marker

p1 <- FeaturePlot(mixed_2, features = B.marker)
p2 <- DotPlot(mixed_2, features = B.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = B.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到2.4 B细胞marker

p1 <- FeaturePlot(mixed_2, features = T.marker)
p2 <- DotPlot(mixed_2, features = T.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到2.5 T细胞marker

p1 <- FeaturePlot(mixed_2, features = DC.marker)
p2 <- DotPlot(mixed_2, features = DC.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = DC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到2.6 DC细胞marker

p1 <- FeaturePlot(mixed_2, features = Granulocyte.marker)
p2 <- DotPlot(mixed_2, features = Granulocyte.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = Granulocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
##得到2.7 粒细胞marker



####细胞注释Mixed（cluster1）##########################
#handcelltype
subtype <- c("0" = "DC",
             "1" = "DC",
             "2" = "B Lymphocyte",
             "3" = "Macrophage")
names(subtype) <- levels(mixed_2)
mixed_2 <- RenameIdents(mixed_2, subtype)
mixed_2@meta.data$subtype <- Idents(mixed_2)
Idents(mixed_2) <- mixed_2$subtype


MVF.marker <- unique(c(Macrophage.marker, B.marker, T.marker, DC.marker)) 
p1 <- FeaturePlot(mixed_2, features = MVF.marker, raster=FALSE)
p2 <- DotPlot(mixed_2, features = MVF.marker) + coord_flip()
p3 <- VlnPlot(mixed_2, features = MVF.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到2.8 Mixed（cluster1）的marker

Idents(mouse_all, cells = colnames(mixed_2)) <- Idents(mixed_2)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()



####筛选Mixed（cluster2）##########################
#筛选Mixed（cluster2）
mixed_3 <- subset(mouse_all, handcelltype == "Erythrocyte")
object_for_clustree <- mixed_3
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到3.1 Mixed（cluster2）分群树状图

mixed_3 <- FindClusters(mixed_3, resolution = 1.0)
mixed_3 <- RunUMAP(mixed_3, dims = 1:25)

DimPlot(mixed_3,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到3.2 UMAP-Mixed（cluster2）分群



####cluster2的marker可视化##########################
#cluster2的marker可视化
p1 <- FeaturePlot(mixed_3, features = Erythrocyte.marker)
p2 <- DotPlot(mixed_3, features = Erythrocyte.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = Erythrocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到3.3 红细胞marker

p1 <- FeaturePlot(mixed_3, features = T.marker)
p2 <- DotPlot(mixed_3, features = T.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到3.4 T细胞marker



####细胞注释Mixed（cluster2）##########################
#handcelltype
subtype <- c("0" = "Erythrocyte",
             "1" = "Erythrocyte")
names(subtype) <- levels(mixed_3)
mixed_3 <- RenameIdents(mixed_3, subtype)
mixed_3@meta.data$subtype <- Idents(mixed_3)
Idents(mixed_3) <- mixed_3$subtype

Idents(mouse_all, cells = colnames(mixed_3)) <- Idents(mixed_3)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()



####筛选cluster3##########################
#筛选cluster3
mixed_4 <- subset(mouse_all, handcelltype == "Mixed（Erythrocyte/B）")
object_for_clustree <- mixed_4
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到4.1 Mixed（cluster3）分群树状图

mixed_4 <- FindClusters(mixed_4, resolution = 0.8)
mixed_4 <- RunUMAP(mixed_4, dims = 1:25)

DimPlot(mixed_4,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到4.2 UMAP-Mixed（cluster3）分群



####cluster3的marker可视化##########################
#cluster3的marker可视化
p1 <- FeaturePlot(mixed_4, features = Erythrocyte.marker)
p2 <- DotPlot(mixed_4, features = Erythrocyte.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = Erythrocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到4.3 红细胞细胞marker

p1 <- FeaturePlot(mixed_4, features = B.marker)
p2 <- DotPlot(mixed_4, features = B.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = B.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到4.4 B细胞marker



####细胞注释Mixed（cluster3）##########################
#handcelltype
subtype <- c("0" = "B Lymphocyte",
             "1" = "Erythrocyte")
names(subtype) <- levels(mixed_4)
mixed_4 <- RenameIdents(mixed_4, subtype)
mixed_4@meta.data$subtype <- Idents(mixed_4)
Idents(mixed_4) <- mixed_4$subtype


p1 <- FeaturePlot(mixed_4, features = c(Erythrocyte.marker, B.marker),
                  raster=FALSE)
p2 <- DotPlot(mixed_4, features = c(Erythrocyte.marker, B.marker))+
  coord_flip()
p3 <- VlnPlot(mixed_4, features = c(Erythrocyte.marker, B.marker),
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到4.5 Mixed（cluster3）的marker

Idents(mouse_all, cells = colnames(mixed_4)) <- Idents(mixed_4)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()



####筛选cluster4##########################
#筛选cluster4
mixed_5 <- subset(mouse_all, handcelltype == "Mixed（MAC/G）")
object_for_clustree <- mixed_5
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到5.1 Mixed（cluster4）分群树状图

mixed_5 <- FindClusters(mixed_5, resolution = 1.0)
mixed_5 <- RunUMAP(mixed_5, dims = 1:25)

DimPlot(mixed_5,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到5.2 UMAP-Mixed（cluster4）分群



####cluster4的marker可视化##########################
#cluster4的marker可视化
p1 <- FeaturePlot(mixed_5, features = Macrophage.marker)
p2 <- DotPlot(mixed_5, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(mixed_5, features = Macrophage.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到5.3 巨噬细胞细胞marker

p1 <- FeaturePlot(mixed_5, features = Granulocyte.marker)
p2 <- DotPlot(mixed_5, features = Granulocyte.marker)+coord_flip()
p3 <- VlnPlot(mixed_5, features = Granulocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到5.4 粒细胞marker



####细胞注释Mixed（cluster4）##########################
#handcelltype
subtype <- c("0" = "Granulocyte",
             "1" = "Macrophage")
names(subtype) <- levels(mixed_5)
mixed_5 <- RenameIdents(mixed_5, subtype)
mixed_5@meta.data$subtype <- Idents(mixed_5)
Idents(mixed_5) <- mixed_5$subtype


p1 <- FeaturePlot(mixed_5, features = c(Granulocyte.marker, Macrophage.marker),
                  raster=FALSE)
p2 <- DotPlot(mixed_5, features = c(Granulocyte.marker, Macrophage.marker))+
  coord_flip()
p3 <- VlnPlot(mixed_5, features = c(Granulocyte.marker, Macrophage.marker),
              pt.size = 0, layer = "counts",
              log = TRUE )
p1 / p2 / p3
  ##得到5.5 Mixed（cluster4）的marker



####免疫细胞亚群写入###############################
#免疫细胞亚群写入
Idents(mouse_all, cells = colnames(mixed_5)) <- Idents(mixed_5)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()
  ##得到6.免疫细胞亚群注释

saveRDS(mouse_all, "./data source/rds/免疫细胞（注释）.rds")

mouse_all <- readRDS("./data source/rds/确定注释.rds")
sub_immune <- readRDS("./data source/rds/免疫细胞（注释）.rds")
Idents(sub_immune)

Idents(mouse_all, cells = colnames(sub_immune)) <- Idents(sub_immune)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5) +
  NoLegend()



####设定细胞顺序###############
#设定细胞顺序
mouse_all$usesubtype <- Idents(mouse_all)
mouse_all$usesubtype <- factor(x = mouse_all$usesubtype,
                            levels = c("EC", "Macrophage", "DC",
                                       "B Lymphocyte", "Granulocyte", 
                                       "Erythrocyte", "PTC", "POD", 
                                       "MC", "FIB", "VSMC"))
levels(mouse_all) <- c("EC", "Macrophage", "DC", "B Lymphocyte",
                       "Granulocyte", "Erythrocyte", "PTC", "POD", 
                       "MC", "FIB", "VSMC")
saveRDS(mouse_all, "./data source/rds/确定注释（免疫）.rds")
