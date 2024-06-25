rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(monocle))BiocManager::install("monocle")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(tidydr))install.packages("tidydr")
  if(!require(ggforce))install.packages("ggforce")
  if(!require(ggrastr))install.packages("ggrastr")
  if(!require(viridis))install.packages("viridis")
  if(!require(ggridges))install.packages("ggridges")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(data.table))install.packages("data.table")
}

####创建monocle文件########################
#创建monocle文件
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_monocle <- readRDS("./data source/rds/POD时序（起点）.rds")



####分支图#################
#分支图
plot_cell_trajectory(mouse_monocle, cell_size = 2.2, color_by = "group") +
  facet_wrap(~group, nrow = 2)
  ##得到1.1 分支图（各分组）

plot_cell_trajectory(mouse_monocle, cell_size = 2.2, color_by = "group")
  ##得到1.2 分支图（分组合并）

plot_cell_trajectory(mouse_monocle, color_by = "Pseudotime")
  ##得到1.3 分支图（时序）

plot_cell_trajectory(mouse_monocle, cell_size = 2.2, color_by = "State") +
  facet_wrap(~State, nrow = 3)
  ##得到1.4 分支图（各时序阶段）

plot_cell_trajectory(mouse_monocle, color_by = "State")
  ##得到1.5 分支图（时序阶段合并）

plot_cell_trajectory(mouse_monocle, color_by = "orig.ident", cell_link_size = 1.5)
  ##得到1.6 分支图（来源）


##提取数据
data_df <- t(reducedDimS(mouse_monocle)) %>% as.data.frame() %>% 
  select_(Component_1 = 1, Component_2 = 2) %>% 
  rownames_to_column("cells") %>% 
  mutate(pData(mouse_monocle)$State) %>% 
  mutate(pData(mouse_monocle)$Pseudotime, 
         pData(mouse_monocle)$orig.ident, 
         pData(mouse_monocle)$group)
colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","group")

##提取轨迹
dp_mst <- minSpanningTree(mouse_monocle)
reduced_dim_coords <- reducedDimK(mouse_monocle)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))
edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% 
              select_(source = "sample_name", 
                      source_prin_graph_dim_1 = "prin_graph_dim_1", 
                      source_prin_graph_dim_2 = "prin_graph_dim_2"), 
            by = "source") %>% 
  left_join(ica_space_df %>% 
              select_(target = "sample_name", 
                      target_prin_graph_dim_1 = "prin_graph_dim_1", 
                      target_prin_graph_dim_2 = "prin_graph_dim_2"), 
            by = "target")


##计算细胞比例
Cellratio <- prop.table(table(data_df$State, data_df$group), 
                        margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"group","Freq")


##做图
gg1 <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =Pseudotime)) + 
  scale_color_viridis()+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_arc(arrow = arrow(length = unit(0.1, "inches"), 
                         type = "closed", angle=30),
           aes(x0 = 0, y0 = 3, r = 3,  start = 0.4*pi,  end = -0.4*pi), lwd = 0.6)+
  geom_arc_bar(data = subset(Cellratio, State == '1'), stat = "pie",
               aes(x0 = 6, y0 = 0, r0 = 0, r = 1.5, 
                   amount = Freq, fill = group)) +
  annotate("text", x = 6, y = 2.5, label = "State 1", size = 4, hjust = 0.5) +
  geom_arc_bar(data = subset(Cellratio, State =='2'),stat = "pie",
               aes(x0 = -5, y0 = 9, r0 = 0, r = 1.5,
                   amount = Freq, fill = group)) +
  annotate("text", x = -5, y = 6.5, label = "State 2", size = 4, hjust = 0.5) +
  geom_arc_bar(data = subset(Cellratio, State=='3'),stat = "pie",
               aes(x0 = -5, y0 = -8, r0 = 0, r = 1.5,
                   amount = Freq, fill = group)) +
  annotate("text", x = -5, y = -5.5, label = "State 3", size = 4, hjust = 0.5) +
  scale_fill_manual(values = c('#67A9CC', '#DA8A87'))
gg1
  ##得到1.7 分支图（分组-时序）


plot_cell_trajectory(mouse_monocle, markers = "Klk8",
                     use_color_gradient=T, cell_size = 1, cell_link_size = 1.5)
  ##得到1.8 分支图（KLK8）


genes <- c("Klk1", "Klk14", "Klk1b22", "Klk1b24", "Klk1b3", "Klk1b4", 
           "Klk1b5", "Klk1b8", "Klk7", "Klk8", "Lif", "Lifr")
genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(mouse_monocle)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
pData(mouse_monocle) <- cbind(pData(mouse_monocle), gene_exp)

data <- pData(mouse_monocle)
colnames(data)
data <- data %>% 
  select("group", "Pseudotime", "Klk1", "Klk14", "Klk1b22", "Klk1b24", 
         "Klk1b3", "Klk1b4", "Klk1b5", "Klk1b8", "Klk7", "Klk8", "Lif", "Lifr")

data_long_m <- melt(data, id.vars = c("group", "Pseudotime"), 
                    measure.vars = 3:14,
                    variable.name = c('gene'),
                    value.name = 'value')
colnames(data_long_m)

ggplot(data_long_m, aes(x = Pseudotime, y = value, color = group))+
  geom_smooth(aes(fill = group))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~gene, scales = "free_y", ncol = 4)+ 
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14))+ 
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))
  ##得到1.9 基因表达时序图



####查看随着pseudotime变化的基因#############################################################################
#查看随着pseudotime变化的基因
pseudotimegenes <- differentialGeneTest(mouse_monocle,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(pseudotimegenes, "随Pseudotime变化的基因.csv")

pseudotimegenes_sig <- pseudotimegenes %>% 
  filter(pval < 0.05) %>%  
  top_n(400, num_cells_expressed) 
Time_genes <- pseudotimegenes_sig %>% 
  pull(gene_short_name) %>% 
  as.character()
p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T)
p
  ##得到2.1 时序模块热图
hp.genes <- as.data.frame(p$tree_row$labels[p$tree_row$order])
write.csv(hp.genes, "时序模块基因.csv")


plotdf <- pData(mouse_monocle)
  ##pData就是seurat可用的mate.data
plotdf$group <- factor(plotdf$group, 
                       levels = c("CON", "DKD"))

p5 <- ggplot(plotdf, aes(x = Pseudotime,
                         y = group,
                         fill = group))+
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8))+
  scale_x_continuous(position = 'top')
p5
  ##得到2.2 时序模块山脊图
