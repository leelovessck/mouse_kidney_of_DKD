rm(list = ls())
gc()

#载入必要的package
{
if(!require(clusterProfiler))install.packages("clusterProfiler")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(viridis))install.packages("viridis")
if(!require(ggraph))install.packages("ggraph")
if(!require(ggsignif))install.packages("ggsignif")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
gc()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")
Idents(mouse_all) <- mouse_all$usetype



####比例条形图（group）#################
#比例条形图（group）
Idents(mouse_all) <- mouse_all$usetype
cellnum <- as.data.frame(table(mouse_all$usetype))
  

mouse_all$group.celltype <- paste(mouse_all$usetype,
                                  mouse_all$group, 
                                  sep = "_")
cellnum.group <- as.data.frame(table(mouse_all$group.celltype))

cellratio_group <- prop.table(table(Idents(mouse_all),mouse_all$group),
                              margin = 2)
cellratio_group <- as.data.frame(cellratio_group)

colourCount = length(unique(cellratio_group$Var1))
p1 <- ggplot(cellratio_group) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),
           stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", 
                             "#882E72", "#FF7F00", "#E78AC3", "#33A02C", 
                             "#B2DF8A", "#A6761D", "#999999", "#1e90ff"))+
  theme_classic() +
  labs(x='Group',y = 'Cell Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    size=0.5, linetype="solid"))
  ##得到1.不同类型细胞比例图（group）



####比例条形图（group-详细）#################
#比例条形图（group-详细）
Idents(mouse_all) <- mouse_all$usesubtype
cellnum <- as.data.frame(table(mouse_all$usesubtype))

mouse_all$group.celltype <- paste(mouse_all$usesubtype,
                                  mouse_all$group, 
                                  sep = "_")
cellnum.group <- as.data.frame(table(mouse_all$group.celltype))

cellratio_group <- prop.table(table(Idents(mouse_all),mouse_all$group),
                              margin = 2)
cellratio_group <- as.data.frame(cellratio_group)

colourCount = length(unique(cellratio_group$Var1))
p1 <- ggplot(cellratio_group) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),
           stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", 
                             "#882E72", "#FF7F00", "#E78AC3", "#33A02C", 
                             "#B2DF8A", "#A6761D", "#999999", "#1e90ff"))+
  theme_classic() +
  labs(x='Group',y = 'Cell Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    size=0.5, linetype="solid"))
  ##得到2.不同类型细胞比例图（group-详细）
