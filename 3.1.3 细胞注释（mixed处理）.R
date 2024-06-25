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
mouse_all <- readRDS("./data source/rds/注释（含mixed）.rds")
Idents(mouse_all) <- mouse_all$handcelltype



####筛选Mixed（EC/Immune Cell）##########################
#筛选Mixed（EC/Immune Cell）
mixed_1 <- subset(mouse_all, handcelltype == "Mixed（EC/Immune Cell）")
object_for_clustree <- mixed_1
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到1.1 Mixed（EC Immune Cell）分群树状图

mixed_1 <- FindClusters(mixed_1, resolution = 0.5)
mixed_1 <- RunUMAP(mixed_1, dims = 1:25)

DimPlot(mixed_1,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到1.2 UMAP-Mixed（EC Immune Cell）分群



####EC/Immune Cell的marker可视化##########################
#EC/Immune Cell的marker可视化
EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mixed_1, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_1, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.3 内皮细胞marker（MAIN）

EC.marker <- c("Egf17", "Cdh5", "Ptprb", "Flt1", "Pecam1", "Cd300lg", 
               "Cyyr1", "Rasgrp3", "F8", "Gpihbp1", "Kdr", "Emcn", 
               "Ehd3", "Srgn", "Egfl7", "S100a10", "Clic4", "Ly6a")
p1 <- FeaturePlot(mixed_1, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_1, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.4 内皮细胞marker（EXTENDED）


IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_1, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_1, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.5 免疫细胞marker（MAIN）

IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_1, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_1, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_1, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.6 免疫细胞marker（EXTENDED）



####细胞注释Mixed（EC/Immune Cell）##########################
#handcelltype
subtype <- c("0" = "EC",
             "1" = "EC",
             "2" = "Immune Cell")
names(subtype) <- levels(mixed_1)
mixed_1 <- RenameIdents(mixed_1, subtype)
mixed_1@meta.data$subtype <- Idents(mixed_1)
Idents(mixed_1) <- mixed_1$subtype

EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_1, features = c(IMMUNE.marker, EC.marker),raster=FALSE)
p2 <- DotPlot(mixed_1, features = c(IMMUNE.marker, EC.marker))+coord_flip()
p3 <- VlnPlot(mixed_1, features = c(IMMUNE.marker, EC.marker),
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到1.7 Mixed（EC Immune Cell）的marker

Idents(mouse_all, cells = colnames(mixed_1)) <- Idents(mixed_1)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()



####筛选Mixed（MC/FIB/VSMC）##########################
#筛选Mixed（MC/FIB/VSMC）
mixed_2 <- subset(mouse_all, handcelltype == "Mixed（MC/FIB/VSMC）")
object_for_clustree <- mixed_2
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到2.1 Mixed（MC FIB VSMC）分群树状图

mixed_2 <- FindClusters(mixed_2, resolution = 0.2)
mixed_2 <- RunUMAP(mixed_2, dims = 1:25)

DimPlot(mixed_2,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到2.2 UMAP-Mixed（MC FIB VSMC）分群



####MC/FIB/VSMC的marker可视化##########################
#MC/FIB/VSMC的marker可视化
MC.marker <- c("Gata3", "Pdgfrb", "Myl9", "Acta2", "Ren1")
p1 <- FeaturePlot(mixed_2, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mixed_2, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.3 系膜细胞marker（MAIN）

MC.marker <- c("Gucy1a3", "Gata3", "Pdgfrb", "Sfrp2", "Itga8", "Myl9", 
               "Agtr1a", "Hopx", "Ptn", "Mgp", "Acta2", "Ren1", "Tpm2",
               "Cald1", "Nr4a2", "Actn1", "Akap12", "Rock1", "Sh3bgrl")
p1 <- FeaturePlot(mixed_2, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mixed_2, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.4 系膜细胞marker（EXTENDED）


FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
p1 <- FeaturePlot(mixed_2, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(mixed_2, features = FIB.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.5 成纤维细胞marker


VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")
p1 <- FeaturePlot(mixed_2, features = VSMC.marker,raster=FALSE)
p2 <- DotPlot(mixed_2, features = VSMC.marker)+coord_flip()
p3 <- VlnPlot(mixed_2, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.6 血管平滑肌细胞marker



####细胞注释Mixed（MC/FIB/VSMC）##########################
#handcelltype
subtype <- c("0" = "MC",
             "1" = "FIB",
             "2" = "VSMC",
             "3" = "MC")
names(subtype) <- levels(mixed_2)
mixed_2 <- RenameIdents(mixed_2, subtype)
mixed_2@meta.data$subtype <- Idents(mixed_2)
Idents(mixed_2) <- mixed_2$subtype


MVF.marker <- unique(c(MC.marker, VSMC.marker, FIB.marker)) 
p1 <- FeaturePlot(mixed_2, features = MVF.marker, raster=FALSE)
p2 <- DotPlot(mixed_2, features = MVF.marker) + coord_flip()
p3 <- VlnPlot(mixed_2, features = MVF.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到2.7 Mixed（MC FIB VSMC）的marker

Idents(mouse_all, cells = colnames(mixed_2)) <- Idents(mixed_2)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()



####筛选Mixed（Unknown）##########################
#筛选Mixed（Unknown）
mixed_3 <- subset(mouse_all, handcelltype == "Mixed（Unknown）")
object_for_clustree <- mixed_3
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到3.1 Mixed（Unknown）分群树状图

mixed_3 <- FindClusters(mixed_3, resolution = 0.6)
mixed_3 <- RunUMAP(mixed_3, dims = 1:25)

DimPlot(mixed_3,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到3.2 UMAP-Mixed（Unknown）分群



####Unknown的marker可视化##########################
#Unknown的marker可视化
EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mixed_3, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.3 内皮细胞marker（MAIN）

EC.marker <- c("Egf17", "Cdh5", "Ptprb", "Flt1", "Pecam1", "Cd300lg", 
               "Cyyr1", "Rasgrp3", "F8", "Gpihbp1", "Kdr", "Emcn", 
               "Ehd3", "Srgn", "Egfl7", "S100a10", "Clic4", "Ly6a")
p1 <- FeaturePlot(mixed_3, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.4 内皮细胞marker（EXTENDED）


MC.marker <- c("Gata3", "Pdgfrb", "Myl9", "Acta2", "Ren1")
p1 <- FeaturePlot(mixed_3, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.5 系膜细胞marker（MAIN）

MC.marker <- c("Gucy1a3", "Gata3", "Pdgfrb", "Sfrp2", "Itga8", "Myl9", 
               "Agtr1a", "Hopx", "Ptn", "Mgp", "Acta2", "Ren1", "Tpm2",
               "Cald1", "Nr4a2", "Actn1", "Akap12", "Rock1", "Sh3bgrl")
p1 <- FeaturePlot(mixed_3, features = MC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = MC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.6 系膜细胞marker（EXTENDED）


POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
p1 <- FeaturePlot(mixed_3, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = POD.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.7 足细胞marker（MAIN）

POD.marker <- c("Nphs1", "Enpep", "Thsd7a", "Dpp4", "Npnt", "Plce1", 
                "Ptpro", "Nphs2", "Clic3", "Cdkn1c", "Synpo", "Mafb", 
                "Wt1", "Cd2ap", "Podxl", "Col4a3", "Golim4", "Magi2", 
                "Pard3b")
p1 <- FeaturePlot(mixed_3, features = POD.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = POD.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.8 足细胞marker（EXTENDED）


IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_3, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.9 免疫细胞marker（MAIN）

IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_3, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.10 免疫细胞marker（EXTENDED）


PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
p1 <- FeaturePlot(mixed_3, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.11 壁上皮细胞marker（MAIN）


PTC.marker <- c("Tacstd2", "Pcbd1", "Cdh16", "Pdzk1ip1", "Atp1b1", "Spp1", 
                "Fxyd2", "Cryab", "Ldhb", "S100g", "Slc12a3")
p1 <- FeaturePlot(mixed_3, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.12 壁上皮细胞marker（EXTENDED）


FIB.marker <- c("Clca3a1", "Dkk2", "Itgbl1", "Col8a1", "Acta2", "Adamtsl1")
p1 <- FeaturePlot(mixed_3, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = FIB.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.13 成纤维细胞marker


VSMC.marker <- c("Acta2", "Myh11", "Cnn1", "Cnn2", "Tagln", "Des", "Smtn")
p1 <- FeaturePlot(mixed_3, features = VSMC.marker,raster=FALSE)
p2 <- DotPlot(mixed_3, features = VSMC.marker)+coord_flip()
p3 <- VlnPlot(mixed_3, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.14 血管平滑肌细胞marker



####细胞注释Mixed（Unknown）##########################
#handcelltype
subtype <- c("0" = "EC/IMMUNE",
             "1" = "PTC",
             "2" = "POD")
names(subtype) <- levels(mixed_3)
mixed_3 <- RenameIdents(mixed_3, subtype)
mixed_3@meta.data$subtype <- Idents(mixed_3)
Idents(mixed_3) <- mixed_3$subtype


EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
MERGE.marker <- c(EC.marker, IMMUNE.marker, PTC.marker, POD.marker)
p1 <- FeaturePlot(mixed_3, features = MERGE.marker, raster=FALSE)
p2 <- DotPlot(mixed_3, features = MERGE.marker) + coord_flip()
p3 <- VlnPlot(mixed_3, features = MERGE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.15 Mixed（Unknown）的marker



####细胞注释Mixed（Unknown亚群）##########################
#handcelltype
mixed_4 <- subset(mixed_3, subtype == "EC/IMMUNE")
object_for_clustree <- mixed_4
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到3.16 Mixed（Unknown亚群）分群树状图

mixed_4 <- FindClusters(mixed_4, resolution = 0.5)
mixed_4 <- RunUMAP(mixed_4, dims = 1:25)

DimPlot(mixed_4,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)

EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
p1 <- FeaturePlot(mixed_4, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_4, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.17 内皮细胞marker（MAIN）

EC.marker <- c("Egf17", "Cdh5", "Ptprb", "Flt1", "Pecam1", "Cd300lg", 
               "Cyyr1", "Rasgrp3", "F8", "Gpihbp1", "Kdr", "Emcn", 
               "Ehd3", "Srgn", "Egfl7", "S100a10", "Clic4", "Ly6a")
p1 <- FeaturePlot(mixed_4, features = EC.marker,raster=FALSE)
p2 <- DotPlot(mixed_4, features = EC.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.18 内皮细胞marker（EXTENDED）


IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_4, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_4, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.19 免疫细胞marker（MAIN）

IMMUNE.marker <- c("Mpeg1", "Cybb", "Lilrb4a", "Ptprc", "Ctss", "Fcer1g", 
                   "Tyrobp", "Lyz2", "Cd74", "H2-Aa")
p1 <- FeaturePlot(mixed_4, features = IMMUNE.marker,raster=FALSE)
p2 <- DotPlot(mixed_4, features = IMMUNE.marker)+coord_flip()
p3 <- VlnPlot(mixed_4, features = IMMUNE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.20 免疫细胞marker（EXTENDED）


newtype <- c("0" = "Immune Cell",
             "1" = "EC")
names(newtype) <- levels(mixed_4)
mixed_4 <- RenameIdents(mixed_4, newtype)
mixed_4@meta.data$newtype <- Idents(mixed_4)
Idents(mixed_4) <- mixed_4$newtype
Idents(mixed_3, cells = colnames(mixed_4)) <- Idents(mixed_4)
DimPlot(mixed_3, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()


####整合Mixed（Unknown）##########################
#整合Mixed（Unknown）
EC.marker <- c("Cdh5", "Flt1", "Pecam1", "Emcn", "Ehd3")
IMMUNE.marker <- c("Lyz2", "Cd74", "H2-Aa")
PTC.marker <- c("Atp1b1", "Slc12a3", "Fxyd2")
POD.marker <- c("Nphs1", "Thsd7a", "Nphs2", "Synpo", "Mafb", "Wt1",
                "Cd2ap", "Podxl", "Col4a3")
MERGE.marker <- c(EC.marker, IMMUNE.marker, PTC.marker, POD.marker)
p1 <- FeaturePlot(mixed_3, features = MERGE.marker, raster=FALSE)
p2 <- DotPlot(mixed_3, features = MERGE.marker) + coord_flip()
p3 <- VlnPlot(mixed_3, features = MERGE.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p1 / p2 / p3
  ##得到3.21 Mixed（Unknown）的marker（新）


Idents(mouse_all, cells = colnames(mixed_3)) <- Idents(mixed_3)
DimPlot(mouse_all, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()



####设定细胞顺序###############
#设定细胞顺序
mouse_all$usetype <- Idents(mouse_all)
mouse_all$usetype <- factor(x = mouse_all$usetype,
                            levels = c("EC", "Immune Cell", "PTC",
                                       "POD", "MC", "FIB", "VSMC"))
levels(mouse_all) <- c("EC", "Immune Cell", "PTC", "POD", "MC", "FIB", "VSMC")


saveRDS(mouse_all, "./data source/rds/确定注释.rds")
