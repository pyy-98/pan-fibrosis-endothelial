library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(circlize)

Heart_function <- read.table("heart_function_celltype.txt" , sep = "\t" , header = T)
Kidney_function <- read.table("kidney_function_celltype.txt" , sep = "\t" , header = T)
Lung_function <- read.table("lung_function_celltype.txt" , sep = "\t" , header = T)
result <- rbind(Lung_function , Heart_function , Kidney_function)

Heart_log10p <- read.table("heart_function_celltype_log10p.txt" , sep = "\t" , header = T)
Kidney_log10p <- read.table("kidney_function_celltype_log10p.txt" , sep = "\t" , header = T)
Lung_log10p <- read.table("lung_function_celltype_log10p.txt" , sep = "\t" , header = T)
result_log10p <- rbind(Lung_log10p , Heart_log10p , Kidney_log10p)

result_log10p[result_log10p < -log10(0.05)] <- ""
result_log10p[result_log10p >= -log10(0.05) & result_log10p < -log10(0.01)] <- "*"
result_log10p[result_log10p >= -log10(0.01)] <- "**"


range(result)
# 设置热图颜色
col_fun <- colorRamp2(c(-0.3 , 0 , 0.3) , c("#6799A5" , "white" , "#C1687B"))
round_cell <-function(j,i,x,y,w,h,fill){
  grid.roundrect(x,y,w,h,r = unit(0.3,'snpc'),
                 gp = gpar(col ="white",fill = fill,lwd =2))
  grid.text(result_log10p[j,i] , x , y , vjust = 0.7 , 
            gp = gpar(fontsize = 13 , col = "white"))
}
# 创建基因分组
gene_groups <- c(rep(1:3 , c(3 , 2 , 2)))
# 自定义颜色
group_colors <- c("#E5D2A2","#DAA1AA","#DCA68F")

# 创建顶部注释（基因分组）
top_annotation <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = "white" , col = "white") ,
    height = unit(0.25 , "cm") , 
    labels = c("Lung" , "Heart" , "Kidney") , 
    labels_gp = gpar(col = "black" , fontsize = 7)
  ) , 
  Group = anno_block(
    gp = gpar(fill = group_colors , col = "white") , 
    height = unit(0.1 , "cm") , 
    labels = NULL , 
    labels_gp = gpar(col = "black" , fontsize = 6)
  ) , 
  show_annotation_name = FALSE
)
result <- t(result)

Heatmap(
  result , 
  name = "Expression" , 
  col = col_fun , 
  # 设置格子样式
  cell_fun = round_cell , 
  rect_gp = gpar(type = "none") , 
  # 设置格子尺寸
  width = ncol(result) * unit(0.6 , "cm") , 
  height = nrow(result) * unit(0.6 , "cm") , 
  # 热图外观设置
  row_names_side = "left" , 
  row_names_gp = gpar(fontsize = 7) , 
  column_names_gp = gpar(fontsize = 7) , 
  column_names_rot = 45 ,
  
  # 使用圆角矩形（通过 raster 参数实现）
  raster_quality = 2 , 
  # 聚类设置（可选）
  cluster_rows = F , 
  cluster_columns = F , 
  # 列分割（用于分组显示）
  column_split = gene_groups , 
  column_title = NULL , 
  
  # 顶部注释
  top_annotation = top_annotation , 
  
  # 颜色条设置
  heatmap_legend_param = list(
    title = "Pearson" , 
    title_gp = gpar(fontsize = 9) , 
    labels_gp = gpar(fontsize = 7)
  )
)

