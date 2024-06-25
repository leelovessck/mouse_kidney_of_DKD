rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
}

setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
getwd()
mouse_all <- readRDS("./data source/rds/确定注释（免疫）.rds")



####EC############################
#EC
EC <- subset(mouse_all, usetype == "EC")
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
table(EC$usetype)
EC_DEG <- FindMarkers(EC, 
                      min.pct = 0.10, 
                      logfc.threshold = 0.10,
                      group.by = "group",
                      ident.1 = "DKD",
                      ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(EC_DEG, "全部EC的DEG.csv")



####PTC############################
#PTC
PTC <- subset(mouse_all, usetype == "PTC")
PTC$usetype <- droplevels(PTC$usetype,
                          exclude = setdiff(
                            levels(PTC$usetype),
                            unique(PTC$usetype)))
table(PTC$usetype)
PTC_DEG <- FindMarkers(PTC, 
                       min.pct = 0.10, 
                       logfc.threshold = 0.10,
                       group.by = "group",
                       ident.1 = "DKD",
                       ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(PTC_DEG, "全部PTC的DEG.csv")



####POD############################
#POD
POD <- subset(mouse_all, usetype == "POD")
POD$usetype <- droplevels(POD$usetype,
                          exclude = setdiff(
                            levels(POD$usetype),
                            unique(POD$usetype)))
table(POD$usetype)
POD_DEG <- FindMarkers(POD, 
                       min.pct = 0.10, 
                       logfc.threshold = 0.10,
                       group.by = "group",
                       ident.1 = "DKD",
                       ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(POD_DEG, "全部POD的DEG.csv")



####MC############################
#MC
MC <- subset(mouse_all, usetype == "MC")
MC$usetype <- droplevels(MC$usetype,
                          exclude = setdiff(
                            levels(MC$usetype),
                            unique(MC$usetype)))
table(MC$usetype)
MC_DEG <- FindMarkers(MC, 
                       min.pct = 0.10, 
                       logfc.threshold = 0.10,
                       group.by = "group",
                       ident.1 = "DKD",
                       ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(MC_DEG, "全部MC的DEG.csv")



####全部细胞##################
#全部细胞
ALL_DEG <- FindMarkers(mouse_all, 
                      min.pct = 0.10, 
                      logfc.threshold = 0.10,
                      group.by = "group",
                      ident.1 = "DKD",
                      ident.2 = "CON")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(ALL_DEG, "全部细胞的DEG.csv")



####各DEG可视化####################
#各DEG可视化
##数据准备
ALL_DEG$celltype <- "All Cell"
EC_DEG$celltype <- "EC"
PTC_DEG$celltype <- "PTC"
POD_DEG$celltype <- "POD"
MC_DEG$celltype <- "MC"

diff_cell<-rbind(ALL_DEG, EC_DEG, PTC_DEG, POD_DEG, MC_DEG)
diff_cell$GENE <- rownames(diff_cell)
head(diff_cell)


##显著性
diff_cell$label <- ifelse(diff_cell$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")

top20_ALL <- filter(diff_cell,celltype=="All Cell") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(20,abs(avg_log2FC))

top20_EC <- filter(diff_cell,celltype=="EC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(20,abs(avg_log2FC))

top20_PTC <- filter(diff_cell,celltype=="PTC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(20,abs(avg_log2FC))

top20_POD <- filter(diff_cell,celltype=="POD") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(20,abs(avg_log2FC))

top20_MC <- filter(diff_cell,celltype=="MC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(20,abs(avg_log2FC))

top20 <- rbind(top20_ALL, top20_EC, top20_PTC, top20_POD, top20_MC)
diff_cell$size <- case_when(!(diff_cell$GENE %in% top20$GENE)~ 1,
                            diff_cell$GENE %in% top20$GENE ~ 2)


##提取非Top10的基因表格；
dt <- filter(diff_cell,size==1)
head(dt)


##绘图
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top20,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

dfbar<-data.frame(x=c("All Cell", "EC", "MC", "POD", "PTC"),
                  y=c(9, 8, 9, 10, 10))
dfbar1<-data.frame(x=c("All Cell", "EC", "MC", "POD", "PTC"),
                   y=c(-8, -7, -9, -10, -10))
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
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top20,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
dfcol<-data.frame(x=c("All Cell", "EC", "MC", "POD", "PTC"),
                  y=0,
                  label=c("All Cell", "EC", "MC", "POD", "PTC"))
p2

mycol <- c("#E64B357F", "#00A0877F", "#34887F", "#F39B7F7F","grey")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
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
    data=top20,
    aes(x=celltype,y=avg_log2FC,label=GENE),
    size =3,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"), 
    max.overlaps = 100
  )
p4
p5 <- p4+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")
p5
light_red <- rgb(250, 100, 100, maxColorValue = 255)
p6 <- p4 +
  scale_color_manual(name=NULL,
                     values = c(light_red,"grey"))+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6
  ##得到9.全部细胞DEG汇总
