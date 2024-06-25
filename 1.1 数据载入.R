rm(list = ls())
gc()

#打开必要的package
{
  if(!require(multtest))BiocManager::install("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(Matrix))install.packages("Matrix")
  if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
  if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(multtest))BiocManager::install("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(clustree))install.packages("clustree")
  if(!require(patchwork))install.packages("patchwork")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()



####载入数据#############################
#载入数据
EXP1 <- read.table("./data source/raw/EXP1.txt",
                   sep = "\t",header = T)
EXP2 <- read.table("./data source/raw/EXP2.txt",
                   sep = "\t",header = T)



####GeneID转化为symbol######################
#GeneID转化为symbol
##EXP1
EXP1$Geneid <- sub("\\..*", "", EXP1$Geneid)  
rownames(EXP1) <- EXP1[,1]
geneid_EXP1 <- rownames(EXP1)
gene_EXP1 <- bitr(geneid_EXP1,
                  fromType = "ENSEMBL",
                  toType = "SYMBOL",
                  OrgDb = org.Mm.eg.db)
colnames(gene_EXP1)[1] <- "Geneid"
symbol_EXP1 <- merge(gene_EXP1, EXP1, by = "Geneid", all = FALSE)  


##EXP2
EXP2$Geneid <- sub("\\..*", "", EXP2$Geneid)  
rownames(EXP2) <- EXP2[,1]
geneid_EXP2 <- rownames(EXP2)
gene_EXP2 <- bitr(geneid_EXP2,
                  fromType = "ENSEMBL",
                  toType = "SYMBOL",
                  OrgDb = org.Mm.eg.db)
colnames(gene_EXP2)[1] <- "Geneid"
symbol_EXP2 <- merge(gene_EXP2, EXP2, by = "Geneid", all = FALSE)  



####数据拆分##########################
#数据拆分
##在每个实验中使用三只对照和三只糖尿病小鼠进行两个独立实验（Exp1和Exp2）
##Exp1包含推荐的Array Control RNA标准品，而Exp 2则省略了该标准品
##Exp1和Exp2的Columns1-10是糖尿病，Columns11-20是对照

##EXP1
colnames(symbol_EXP1)
sample_col_names <- colnames(symbol_EXP1)
s_numbers <- as.numeric(sub(".*S(\\d+)$", "\\1", sample_col_names))  
dkd_cols <- c("SYMBOL", sample_col_names[s_numbers %in% 1:10])  
con_cols <- c("SYMBOL", sample_col_names[s_numbers %in% 11:20])  
EXP1_DKD <- symbol_EXP1[, dkd_cols, drop = FALSE]  
EXP1_CON <- symbol_EXP1[, con_cols, drop = FALSE]  


##EXP2
colnames(symbol_EXP2)
sample_col_names <- colnames(symbol_EXP2)
s_numbers <- as.numeric(sub(".*S(\\d+)$", "\\1", sample_col_names))  
dkd_cols <- c("SYMBOL", sample_col_names[s_numbers %in% 1:10])  
con_cols <- c("SYMBOL", sample_col_names[s_numbers %in% 11:20])  
EXP2_DKD <- symbol_EXP2[, dkd_cols, drop = FALSE]  
EXP2_CON <- symbol_EXP2[, con_cols, drop = FALSE]  



####分别建立seurat对象#######################
#分别建立seurat对象
EXP1_DKD <- EXP1_DKD %>% 
  distinct(SYMBOL, .keep_all = T)
rownames(EXP1_DKD) <- EXP1_DKD[,1]
EXP1_DKD <- EXP1_DKD[,-1]
EXP1_DKD <- as(as.matrix(EXP1_DKD), "dgCMatrix")
EXP1_DKD_Seurat <- CreateSeuratObject(counts = EXP1_DKD,
                                      min.cells = 0,
                                      min.features = 0)
EXP1_DKD_Seurat$group <- "DKD"
EXP1_DKD_Seurat$source <- "EXP1_DKD"
Idents(EXP1_DKD_Seurat)
row.names(EXP1_DKD_Seurat[["RNA"]]) 

EXP1_CON <- EXP1_CON %>% 
  distinct(SYMBOL, .keep_all = T)
rownames(EXP1_CON) <- EXP1_CON[,1]
EXP1_CON <- EXP1_CON[,-1]
EXP1_CON <- as(as.matrix(EXP1_CON), "dgCMatrix")
EXP1_CON_Seurat <- CreateSeuratObject(counts = EXP1_CON)
EXP1_CON_Seurat$group <- "CON"
EXP1_CON_Seurat$source <- "EXP1_CON"
Idents(EXP1_CON_Seurat)
row.names(EXP1_CON_Seurat[["RNA"]]) 

EXP2_DKD <- EXP2_DKD %>% 
  distinct(SYMBOL, .keep_all = T)
rownames(EXP2_DKD) <- EXP2_DKD[,1]
EXP2_DKD <- EXP2_DKD[,-1]
EXP2_DKD <- as(as.matrix(EXP2_DKD), "dgCMatrix")
EXP2_DKD_Seurat <- CreateSeuratObject(counts = EXP2_DKD)
EXP2_DKD_Seurat$group <- "DKD"
EXP2_DKD_Seurat$source <- "EXP2_DKD"
Idents(EXP2_DKD_Seurat)
row.names(EXP2_DKD_Seurat[["RNA"]]) 

EXP2_CON <- EXP2_CON %>% 
  distinct(SYMBOL, .keep_all = T)
rownames(EXP2_CON) <- EXP2_CON[,1]
EXP2_CON <- EXP2_CON[,-1]
EXP2_CON <- as(as.matrix(EXP2_CON), "dgCMatrix")
EXP2_CON_Seurat <- CreateSeuratObject(counts = EXP2_CON)
EXP2_CON_Seurat$group <- "CON"
EXP2_CON_Seurat$source <- "EXP2_CON"
Idents(EXP2_CON_Seurat)
row.names(EXP2_CON_Seurat[["RNA"]]) 



####PCA分析#######################
#PCA分析
EXP1_DKD_Seurat <- NormalizeData(EXP1_DKD_Seurat, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
EXP1_DKD_Seurat <- FindVariableFeatures(EXP1_DKD_Seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
EXP1_DKD_Seurat <- ScaleData(EXP1_DKD_Seurat)
EXP1_DKD_Seurat <- RunPCA(EXP1_DKD_Seurat,
                    features = VariableFeatures(object = EXP1_DKD_Seurat))
EXP1_DKD_Seurat<- JackStraw(EXP1_DKD_Seurat,
                      num.replicate = 100,
                      dims = 40)
EXP1_DKD_Seurat <- ScoreJackStraw(EXP1_DKD_Seurat, dims = 1:40)
p1 <- JackStrawPlot(EXP1_DKD_Seurat, dims = 1:40)
p2 <- ElbowPlot(EXP1_DKD_Seurat, ndims = 50) 
p1 + p2
  ##EXP1_DKD_Seurat使用25个PCA

EXP1_CON_Seurat <- NormalizeData(EXP1_CON_Seurat, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
EXP1_CON_Seurat <- FindVariableFeatures(EXP1_CON_Seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
EXP1_CON_Seurat <- ScaleData(EXP1_CON_Seurat)
EXP1_CON_Seurat <- RunPCA(EXP1_CON_Seurat,
                          features = VariableFeatures(object = EXP1_CON_Seurat))
EXP1_CON_Seurat<- JackStraw(EXP1_CON_Seurat,
                            num.replicate = 100,
                            dims = 40)
EXP1_CON_Seurat <- ScoreJackStraw(EXP1_CON_Seurat, dims = 1:40)
p1 <- JackStrawPlot(EXP1_CON_Seurat, dims = 1:40)
p2 <- ElbowPlot(EXP1_CON_Seurat, ndims = 50) 
p1 + p2
  ##EXP1_CON_Seurat使用25个PCA

EXP2_DKD_Seurat <- NormalizeData(EXP2_DKD_Seurat, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
EXP2_DKD_Seurat <- FindVariableFeatures(EXP2_DKD_Seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
EXP2_DKD_Seurat <- ScaleData(EXP2_DKD_Seurat)
EXP2_DKD_Seurat <- RunPCA(EXP2_DKD_Seurat,
                          features = VariableFeatures(object = EXP2_DKD_Seurat))
EXP2_DKD_Seurat<- JackStraw(EXP2_DKD_Seurat,
                            num.replicate = 100,
                            dims = 40)
EXP2_DKD_Seurat <- ScoreJackStraw(EXP2_DKD_Seurat, dims = 1:40)
p1 <- JackStrawPlot(EXP2_DKD_Seurat, dims = 1:40)
p2 <- ElbowPlot(EXP2_DKD_Seurat, ndims = 50) 
p1 + p2
  ##EXP2_DKD_Seurat使用25个PCA

EXP2_CON_Seurat <- NormalizeData(EXP2_CON_Seurat, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
EXP2_CON_Seurat <- FindVariableFeatures(EXP2_CON_Seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
EXP2_CON_Seurat <- ScaleData(EXP2_CON_Seurat)
EXP2_CON_Seurat <- RunPCA(EXP2_CON_Seurat,
                          features = VariableFeatures(object = EXP2_CON_Seurat))
EXP2_CON_Seurat<- JackStraw(EXP2_CON_Seurat,
                            num.replicate = 100,
                            dims = 40)
EXP2_CON_Seurat <- ScoreJackStraw(EXP2_CON_Seurat, dims = 1:40)
p1 <- JackStrawPlot(EXP2_CON_Seurat, dims = 1:40)
p2 <- ElbowPlot(EXP2_CON_Seurat, ndims = 50) 
p1 + p2
  ##EXP2_CON_Seurat使用25个PCA



####合并seurat对象############################
#合并seurat对象
EXP1_Seurat <- FindIntegrationAnchors(object.list = list(EXP1_CON_Seurat,
                                                         EXP1_DKD_Seurat),
                                      dims = 1:25)
EXP1_integrated <- IntegrateData(anchorset = EXP1_Seurat, dims = 1:25)
row.names(EXP1_integrated[["RNA"]]) 

EXP2_Seurat <- FindIntegrationAnchors(object.list = list(EXP2_CON_Seurat,
                                                         EXP2_DKD_Seurat),
                                      dims = 1:25)
EXP2_integrated <- IntegrateData(anchorset = EXP2_Seurat, dims = 1:25)
row.names(EXP2_integrated[["RNA"]]) 

mouse_anchors <- FindIntegrationAnchors(object.list = list(EXP1_integrated,
                                                           EXP2_integrated),
                                        dims = 1:25)
all_mouse <- IntegrateData(anchorset = mouse_anchors, dims = 1:25)
DefaultAssay(all_mouse) <- "RNA"
row.names(all_mouse[["RNA"]]) 
row.names(all_mouse[["integrated"]]) 



####保存数据#############################
saveRDS(all_mouse, "./data source/rds/质量控制前.rds")
