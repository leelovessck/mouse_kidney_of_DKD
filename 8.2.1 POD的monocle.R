rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(monocle))BiocManager::install("monocle")
}

####创建monocle文件########################
#创建monocle文件
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")

mouse_data <- subset(mouse_all, usetype == "POD")
table(mouse_data$usetype)
mouse_data$usetype <- droplevels(mouse_data$usetype,
                                 exclude = setdiff(
                                   levels(mouse_data$usetype),
                                   unique(mouse_data$usetype)))

mouse_expr <- as.sparse(GetAssayData(mouse_data, layer = "counts"))
mouse_idents <- mouse_data@meta.data
mouse_features <- as.matrix(mouse_expr[,1])
mouse_features[,1] <- row.names(mouse_features)
colnames(mouse_features)[1] <- c("gene_short_name")
mouse_features <- data.frame(mouse_features)

mouse_idents <- new("AnnotatedDataFrame", data = mouse_idents)
mouse_features <- new("AnnotatedDataFrame", data = mouse_features)
mouse_monocle <- newCellDataSet(as.matrix(mouse_expr),
                                phenoData = mouse_idents, 
                                featureData = mouse_features, 
                                expressionFamily=negbinomial.size())

mouse_monocle <- estimateSizeFactors(mouse_monocle)
  ##估计大小因子
mouse_monocle <- estimateDispersions(mouse_monocle)
  ##估计分散性


mouse_monocle <- detectGenes(mouse_monocle, min_expr = 0.05)
  ##至少5%的表达
print(head(fData(mouse_monocle)))
print(head(pData(mouse_monocle)))
expressed_genes <- row.names(subset(fData(mouse_monocle), 
                                    num_cells_expressed >= 5))
  ##至少有5个细胞表达的基因
pData(mouse_monocle)$Total_mRNAs <- Matrix::colSums(exprs(mouse_monocle))
  ##将UMI加入cds



####计算差异基因##################
#计算差异基因
diff_test_res <- differentialGeneTest(mouse_monocle,
                                      fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
mouse_monocle <- setOrderingFilter(mouse_monocle, ordering_genes)
plot_ordering_genes(mouse_monocle)



####降维排序#################
#降维排序
mouse_monocle <- reduceDimension(mouse_monocle, 
                                 max_components = 2,
                                 reduction_method = 'DDRTree')
mouse_monocle <- orderCells(mouse_monocle)
saveRDS(mouse_monocle, "./data source/rds/时序（未定起点）.rds")



####指定起点################
#指定起点
GM_state <- function(cds){
   if (length(unique(pData(cds)$State)) > 1){
     T0_counts <- table(pData(cds)$State, pData(cds)$group)[,"CON"]
     return(as.numeric(names(T0_counts)[which
                                        (T0_counts == max(T0_counts))]))
   } else {
     return (1)
   }
 }
mouse_monocle <- orderCells(mouse_monocle, root_state = GM_state(mouse_monocle))
saveRDS(mouse_monocle, "./data source/rds/POD时序（起点）.rds")



####保存数据######################
#保存数据
plotdf <- pData(mouse_monocle)
mouse_data@meta.data <- plotdf
saveRDS(mouse_data, "./data source/rds/POD时序.rds")
