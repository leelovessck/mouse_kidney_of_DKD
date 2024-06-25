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
object_for_clustree <- sub_immune
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到1.免疫细胞分群树状图

sub_immune <- FindClusters(sub_immune, resolution = 0.5)
sub_immune <- RunUMAP(sub_immune, dims = 1:25)

DimPlot(sub_immune,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
##得到2.UMAP-免疫细胞分群



####细胞marker####################
#细胞marker
markers.all <- FindAllMarkers(sub_immune,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
write.csv(markers.all,file = "免疫细胞各cluster的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "免疫细胞各cluster的markers（top20）.csv")



####GPT调用#################
#GPT调用
res <- gptcelltype(markers.all, 
                   tissuename = 'mouse immune cell', 
                   model = 'gpt-4',
                   topgenenumber = 10)
res
##"Identify cell types of human PBMC cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. \n0:CD7,GNLY,CCL5\n1:S100A8,TYMP,S100A9\n2:HLA-DPB1,MS4A1,HLA-DQB1"

##0: 巨噬细胞（Macrophage）或其他髓系细胞
##1: B细胞（B Lymphocyte）
##2: 红细胞（Erythrocyte）或红细胞前体细胞，这些标记与红细胞相关，但并非免疫细胞特有
##3: B细胞（B Lymphocyte）或浆细胞（Plasma Cell）
##4: 粒细胞（Granulocyte），可能是中性粒细胞（Neutrophil）或其他类型的粒细胞，也可能包含单核细胞（Monocyte）或巨噬细胞（Macrophage）



####细胞注释（AI,TOP20）##########################
#AIcelltype(TOP20)
AIcelltype <- c("0" = "Macrophage",
                "1" = "B Lymphocyte",
                "2" = "Erythrocyte",
                "3" = "B Lymphocyte",
                "4" = "Granulocyte")
names(AIcelltype) <- levels(sub_immune)
sub_immune <- RenameIdents(sub_immune, AIcelltype)
sub_immune@meta.data$AIcelltype <- Idents(sub_immune)
Idents(sub_immune) <- sub_immune$AIcelltype



####可视化（AI,TOP20）##########################
#可视化（AI,TOP20）
DimPlot(sub_immune,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "green", "purple"))
  ##得到2.UMAP-AI注释（TOP20）

table(sub_immune$AIcelltype)
mouse_pearson <- AverageExpression(sub_immune,
                                   group.by = "AIcelltype",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),20000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到3.AI注释的spearman相关性（TOP20）


Macrophage.marker <- c("Stab1", "Apoe", "Lyz2", "Tgfbi", "Vcan")
p1 <- FeaturePlot(sub_immune, features = Macrophage.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Macrophage.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Macrophage.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到4.AI注释巨噬细胞marker（TOP20）

B.marker <- c("Pax5", "Ebf1", "Cd79a", "Cd79b")
p1 <- FeaturePlot(sub_immune, features = B.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = B.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = B.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到5.AI注释B细胞marker（TOP20）

T.marker <- c("Skap1", "Cd247", "Cd3g", "Cd8b1", "Cd8a", "Cd4")
p1 <- FeaturePlot(sub_immune, features = T.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = T.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = T.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到6.AI注释T细胞marker（TOP20）

DC.marker <- c("Clec9a", "Ppt1", "Itgam", "Sirpa", "Itgax", "H2-Aa", 
               "H2-Ab1", "Cd74", "Ccl5")
p1 <- FeaturePlot(sub_immune, features = DC.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = DC.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = DC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到7.AI注释DC细胞marker（TOP20）

Erythrocyte.marker <- c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")
p1 <- FeaturePlot(sub_immune, features = Erythrocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Erythrocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Erythrocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到8.AI注释红细胞marker（TOP20）

Granulocyte.marker <- c("Thbs1", "Cd14", "Junb", "Lcn2", "S100a9", 
                        "Il1b", "Ly6g")
p1 <- FeaturePlot(sub_immune, features = Granulocyte.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = Granulocyte.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = Granulocyte.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到9.AI注释粒细胞marker（TOP20）

NK.marker <- c("Nkg7", "Gzmb", "Gzma")
p1 <- FeaturePlot(sub_immune, features = NK.marker,raster=FALSE)
p2 <- DotPlot(sub_immune, features = NK.marker)+coord_flip()
p3 <- VlnPlot(sub_immune, features = NK.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.AI注释NK细胞marker（TOP20）
