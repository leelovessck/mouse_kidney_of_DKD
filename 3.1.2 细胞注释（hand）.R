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
mouse_all <- readRDS("./data source/rds/样本分群.rds")
Idents(mouse_all) <- mouse_all$seurat_clusters



####marker可视化##########################
#marker可视化
EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mouse_all, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.1 内皮细胞marker（MAIN）

EC.marker <- c("Egf17", "Cdh5", "Ptprb", "Flt1", "Pecam1", "Cd300lg", 
               "Cyyr1", "Rasgrp3", "F8", "Gpihbp1", "Kdr", "Emcn", 
               "Ehd3", "Srgn", "Egfl7", "S100a10", "Clic4", "Ly6a")
p1 <- FeaturePlot(mouse_all, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.2 内皮细胞marker（EXTENDED）


MC.marker <- c("Gata3", "Pdgfrb", "Myl9", "Acta2", "Ren1")
p1 <- FeaturePlot(mouse_all, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.1 系膜细胞marker（MAIN）

MC.marker <- c("Gucy1a3", "Gata3", "Pdgfrb", "Sfrp2", "Itga8", "Myl9", 
               "Agtr1a", "Hopx", "Ptn", "Mgp", "Acta2", "Ren1", "Tpm2",
               "Cald1", "Nr4a2", "Actn1", "Akap12", "Rock1", "Sh3bgrl")
p1 <- FeaturePlot(mouse_all, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.2 系膜细胞marker（EXTENDED）


POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
p1 <- FeaturePlot(mouse_all, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = POD.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.1 足细胞marker（MAIN）

POD.marker <- c("Nphs1", "Enpep", "Thsd7a", "Dpp4", "Npnt", "Plce1", 
                "Ptpro", "Nphs2", "Clic3", "Cdkn1c", "Synpo", "Mafb", 
                "Wt1", "Cd2ap", "Podxl", "Col4a3", "Golim4", "Magi2", 
                "Pard3b")
p1 <- FeaturePlot(mouse_all, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = POD.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.2 足细胞marker（EXTENDED）


IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mouse_all, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到4.1 免疫细胞marker（MAIN）

IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mouse_all, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到4.2 免疫细胞marker（EXTENDED）


PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
p1 <- FeaturePlot(mouse_all, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到5.1 壁上皮细胞marker（MAIN）


PTC.marker <- c("Tacstd2", "Pcbd1", "Cdh16", "Pdzk1ip1", "Atp1b1", "Spp1", 
                "Fxyd2", "Cryab", "Ldhb", "S100g", "Slc12a3")
p1 <- FeaturePlot(mouse_all, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到5.2 壁上皮细胞marker（EXTENDED）


FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
p1 <- FeaturePlot(mouse_all, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = FIB.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到6.成纤维细胞marker


VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")
p1 <- FeaturePlot(mouse_all, features = VSMC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = VSMC.marker)+coord_flip()
p3 <- VlnPlot(mouse_all, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到7.血管平滑肌细胞marker




####细胞注释（hand）##########################
#handcelltype
handcelltype <- c("0" = "EC",
                  "1" = "Mixed（MC/FIB/VSMC）",
                  "2" = "Mixed（Unknown）",
                  "3" = "Mixed（EC/Immune Cell）",
                  "4" = "POD",
                  "5" = "PTC",
                  "6" = "POD",
                  "7" = "Immune Cell",
                  "8" = "EC",
                  "9" = "Immune Cell",
                  "10" = "Mixed（MC/FIB/VSMC）")
names(handcelltype) <- levels(mouse_all)
mouse_all <- RenameIdents(mouse_all, handcelltype)
mouse_all@meta.data$handcelltype <- Idents(mouse_all)
Idents(mouse_all) <- mouse_all$handcelltype
saveRDS(mouse_all, "./data source/rds/注释（含mixed）.rds")



####可视化（hand）##########################
#可视化（hand）
DimPlot(mouse_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE,
        cols = c("red", "blue", "green", "cyan",
                 "yellow", "purple", "orange"))
  ##得到8.UMAP-手动注释

table(mouse_all$handcelltype)
mouse_pearson <- AverageExpression(mouse_all,
                                   group.by = "handcelltype",
                                   assays = "RNA")
mouse_pearson <- mouse_pearson[[1]]
mouse_pearson <- as.matrix(mouse_pearson)
gene_pearson <- names(tail(sort(apply(mouse_pearson, 1, sd)),20000))
cell_pearson <- cor(mouse_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到9.手动注释的spearman相关性

EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mouse_all, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = EC.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.1 手动注释内皮细胞marker

MC.marker <- c("Gata3", "Pdgfrb", "Myl9", "Acta2", "Ren1")
p1 <- FeaturePlot(mouse_all, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = MC.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.2 手动注释系膜细胞marker

POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
p1 <- FeaturePlot(mouse_all, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = POD.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.3 手动注释足细胞marker

IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mouse_all, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = IMMUNE.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.4 手动注释免疫细胞marker

PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
p1 <- FeaturePlot(mouse_all, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = PTC.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.5 手动注释壁上皮细胞marker

VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")
p1 <- FeaturePlot(mouse_all, features = VSMC.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = VSMC.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.6 手动注释血管平滑肌细胞marker

FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
p1 <- FeaturePlot(mouse_all, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(mouse_all, features = FIB.marker)+coord_flip() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- VlnPlot(mouse_all, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到10.7 手动注释成纤维细胞marker
