rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(monocle))BiocManager::install("monocle")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(ggplot2))BiocManager::install("ggplot2")
  if(!require(ggrepel))BiocManager::install("ggrepel")
}

####载入数据###############
#载入数据
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_monocle <- readRDS("./data source/rds/时序（起点）.rds")
monocle_meta.data <- pData(mouse_monocle)

mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")

mouse_data <- subset(mouse_all, usetype == "EC")
table(mouse_data$usetype)
mouse_data$usetype <- droplevels(mouse_data$usetype,
                                 exclude = setdiff(
                                   levels(mouse_data$usetype),
                                   unique(mouse_data$usetype)))

mouse_data@meta.data <- monocle_meta.data

mouse_dkd <- subset(mouse_data, group == "DKD")
table(mouse_dkd$State)

mouse_con <- subset(mouse_data, group == "CON")
table(mouse_con$State)



####计算时序分组marker###################
#计算时序分组marker
Idents(mouse_dkd) <- mouse_dkd$State

markers.dkd <- FindAllMarkers(mouse_dkd,
                              only.pos = TRUE,
                              min.pct = 0.10, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.dkd <- markers.dkd %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
write.csv(markers.dkd, "EC时序分组的marker.csv")


markers.dkd$label <- ifelse(markers.dkd$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")

C1 <- filter(markers.dkd, cluster=="1") %>% 
  distinct(gene,.keep_all = T) %>% 
  top_n(30,abs(avg_log2FC))

C2 <- filter(markers.dkd, cluster=="2") %>% 
  distinct(gene,.keep_all = T) %>% 
  top_n(30,abs(avg_log2FC))

C3 <- filter(markers.dkd, cluster=="3") %>% 
  distinct(gene,.keep_all = T) %>% 
  top_n(30,abs(avg_log2FC))

top <- rbind(C1, C2, C3)
markers.dkd$size <- case_when(!(markers.dkd$gene %in% top$gene)~ 1,
                              markers.dkd$gene %in% top$gene ~ 2)

dt <- markers.dkd %>%  
  filter(size == 1)
head(dt)


##绘图
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

dfbar<-data.frame(x=c("1", "2", "3"),
                  y=c(9, 9, 9))
dfbar1<-data.frame(x=c("1", "2", "3"),
                   y=c(-1, -1, -1))
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
dfcol<-data.frame(x=c("1", "2", "3"),
                  y=0,
                  label=c("1", "2", "3"))
p2

mycol <- c("#E64B357F", "#00A0877F", "grey")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x = x,y = y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3
p4 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)+
  geom_text_repel(
    data=top,
    aes(x=cluster,y=avg_log2FC,label=gene),
    size =3,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"), 
    max.overlaps = 100
  )
p4
p5 <- p4+
  labs(x="cluster",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")
p5
light_red <- rgb(250, 100, 100, maxColorValue = 255)
p6 <- p4 +
  scale_color_manual(name=NULL,
                     values = c(light_red,"grey"))+
  labs(x="cluster",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6
filename <- "3.EC时序分组的marker.png"  
png(filename, width = 1000, height = 750)  
print(p6)  
dev.off()
  ##得到3.EC时序分组的marker



####计算DEG####################
#计算DEG
EC_DEG <- FindMarkers(mouse_dkd, 
                      min.pct = 0.10, 
                      logfc.threshold = 0.10,
                      group.by = "State",
                      ident.1 = "2",
                      ident.2 = "3")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(EC_DEG, "EC时序分组的DEG.csv")


p1 <- EnhancedVolcano(EC_DEG,
                      lab = rownames(EC_DEG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 10e-5,
                      FCcutoff = 1,
                      pointSize = 1.0,
                      labSize = 6.0,
                      title = "EC时序分组的DEG（2倍）")
filename <- "1.EC时序分组的DEG（2倍）.png"  
png(filename, width = 1000, height = 750)  
print(p1)  
dev.off()
  ##得到1.EC时序分组的DEG（2倍）

p1 <- EnhancedVolcano(EC_DEG,
                      lab = rownames(EC_DEG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 10e-5,
                      FCcutoff = 2,
                      pointSize = 1.0,
                      labSize = 6.0,
                      title = "全部EC的DEG（4倍）")
filename <- "2.EC时序分组的DEG（4倍）.png"  
png(filename, width = 1000, height = 750)  
print(p1)  
dev.off()
  ##得到2.EC时序分组的DEG（4倍）
