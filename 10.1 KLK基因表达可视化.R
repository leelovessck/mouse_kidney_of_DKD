rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
  if(!require(ggsignif))install.packages("ggsignif")
  if(!require(cowplot))install.packages("cowplot")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")



####Klk8在各细胞类型的表达#########################
#Klk8在各细胞类型的表达
cell_colors <- c(  
  "Immune Cell" = "#FF0000",
  "EC" = "#3CB44B",          
  "PTC" = "#FFE119",
  "POD" = "#4363D8", 
  "MC" = "#F58231", 
  "FIB" = "#911EBB",     
  "VSMC" = "navy")  
p1 <- DimPlot(mouse_all, reduction = "umap", group.by = "usetype", 
              cols = cell_colors,   
              label = TRUE, pt.size = 0.5) +
  NoLegend()
p2 <- FeaturePlot(mouse_all, features = "Klk8", raster=FALSE)
p3 <- plot_grid(plotlist = list(p1, p2), 
               ncol = 2, 
               labels = c("细胞类型注释", 
                          "Klk8基因的表达"))
p3

p4 <- VlnPlot(mouse_all, 
              features = "Klk8",
              pt.size = 0.1, 
              group.by = "usetype")+
  NoLegend()
p5 <- plot_grid(plotlist = list(p3, p4), 
                ncol = 1)
p5

filenames <- "1.Klk8在各细胞类型的表达.png"
png(filenames, width = 1000, height = 1000)
print(p5)
dev.off()



####EC亚群Klk家族的表达#########################
#EC亚群Klk家族的表达
klk <- c("Klk8", "Klk1", "Klk1b22", "Klk1b24", "Klk1b27", "Klk1b4", 
         "Klk1b5", "Klk1b8", "Klk1b9", "Klk4", "Klk6", "Klk7", "Klk9")
p1 <- DotPlot(mouse_all, features = klk, group.by = "usetype")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

mouse <- read.csv("./result/9.1.1 KLK-基因/长-KLK_EC.csv")
mouse$gene <- factor(mouse$gene, levels = klk)  
p2 <- ggplot(mouse,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene, scales = "fixed", nrow = 1) 
blank_plot <- ggplot() + theme_void() 
p3 <- plot_grid(plotlist = list(p1, blank_plot, p2),   
                ncol = 1,   
                rel_heights = c(1, 0.2, 1))
p3

filenames <- "2. Klk家族的表达.png"
png(filenames, width = 1000, height = 750)
print(p3)
dev.off()
