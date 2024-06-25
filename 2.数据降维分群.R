rm(list = ls())
gc()

#打开必要的package
{
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(multtest))BiocManager::install("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(clustree))install.packages("clustree")
  if(!require(patchwork))install.packages("patchwork")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(biomaRt))BiocManager::install("biomaRt")
}

####载入数据#############################
#载入数据
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/质量控制前.rds")
Idents(mouse_all) <- mouse_all$source
cellcycle <- readRDS("./data source/rds/小鼠细胞周期基因.rds")
mouse_all <- JoinLayers(mouse_all)



####基因质量控制#########################
#基因质量控制
mt.genes <- rownames(mouse_all)[grep("^Mt[sl]",rownames(mouse_all))]
mouse_all <- PercentageFeatureSet(mouse_all, "^Mt[sl]", col.name = "percent.mt")

rb.genes <- rownames(mouse_all)[grep("^Rp[sl]",rownames(mouse_all))]
mouse_all <- PercentageFeatureSet(mouse_all, "^Rp[sl]", col.name = "percent.ribo")


VlnPlot(mouse_all, 
        features = c("nFeature_RNA", "nCount_RNA", 
                     "percent.mt", "percent.ribo"),
        raster=FALSE,
        group.by = "source",
        pt.size = 0,
        ncol = 2)
  ##得到1.基因质量控制图



####细胞周期#########################
#细胞周期
s.genes <- cellcycle$s.genes
g2m.genes <- cellcycle$g2m.genes
mouse_all <- CellCycleScoring(mouse_all, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes,
                              set.ident = TRUE)
VlnPlot(mouse_all,
        features = c("S.Score","G2M.Score"),
        group.by = "source",
        pt.size =0.25)
  ##得到2.细胞周期-小提琴图（source）



####降维###########################
#降维
mouse_all <- ScaleData(mouse_all, features = rownames(mouse_all))
mouse_all <- RunPCA(mouse_all, npcs = 50, verbose = TRUE)
mouse_all<- JackStraw(mouse_all,
                      num.replicate = 100,
                      dims = 40)
mouse_all <- ScoreJackStraw(mouse_all, dims = 1:40)
p1 <- JackStrawPlot(mouse_all, dims = 1:40)
p2 <- ElbowPlot(mouse_all, ndims = 50) 
p1 + p2
  ##得到3.PCA分析图
  ##mouse_all使用25个PCA
mouse_all <- FindNeighbors(mouse_all, dims = 1:25)



####分群############################
#分群
##设置不同resolutions查看效果
object_for_clustree <- mouse_all
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到4.分群树状图


##分群
mouse_all <- FindClusters(mouse_all, resolution = 0.7)
mouse_all <- RunUMAP(mouse_all, dims = 1:25)
saveRDS(mouse_all, "./data source/rds/样本分群.rds")



##细胞marker（cluster）
markers.all <- FindAllMarkers(mouse_all,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
write.csv(markers.all,file = "各cluster的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "各cluster的markers（top20）.csv")



####分群可视化############################
#分群可视化
Idents(mouse_all) <- mouse_all$seurat_clusters
DimPlot(mouse_all,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到5.UMAP-细胞分群

Idents(mouse_all) <- mouse_all$orig.ident
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到6.UMAP-实验次数

Idents(mouse_all) <- mouse_all$group
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到7.UMAP-组别

Idents(mouse_all) <- mouse_all$source
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "orange", "green"))
  ##得到8.UMAP-实验+组别

Idents(mouse_all) <- mouse_all$Phase
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "yellow"))
  ##得到9.UMAP-细胞周期

Idents(mouse_all) <- mouse_all$seurat_clusters
DimPlot(mouse_all,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE, ncol = 2,
        split.by = "source")
  ##得到10.UMAP-细胞分群（source）



####分群的spearman相关########################
#分群的spearman相关
table(mouse_all$seurat_clusters)
mouse_pearson <- AverageExpression(mouse_all,
                                   group.by = "seurat_clusters",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),20000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到11.细胞分群的spearman相关性
