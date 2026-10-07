library(ComplexHeatmap)
library(CellChat)
library(Seurat)
library(ggplot2)
cellchat_Fibrosis <- readRDS("Fibrosis_cellchat.rds")
cellchat_Normal <- readRDS("Normal_cellchat.rds")

object.list <- list(Normal = cellchat_Normal , Fibrosis = cellchat_Fibrosis)
cellchat <- mergeCellChat(object.list, add.names = names(object.list) , cell.prefix = T)

measure <- "count"  ## "weight"
slot.name <- c("netP" , "net")
comparison <- c(1 , 2)
color.heatmap = c("#446d8e", "#993042")
font.size.title = 10
font.size = 8
obj1 <- cellchat@net[[comparison[1]]][[measure]]
obj2 <- cellchat@net[[comparison[2]]][[measure]]
net.diff <- obj2 - obj1
title.name <- NULL
if (measure == "count") {
  if (is.null(title.name)) {
    title.name = "Differential number of interactions"
  }
}else if (measure == "weight") {
  if (is.null(title.name)) {
    title.name = "Differential interaction strength"
  }
}
legend.name = "Relative values"
net <- net.diff
net[is.na(net)] <- 0
mat <- net
if (min(mat) < 0) {
  color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                 c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
  colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                    "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                    na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",max(mat, na.rm = T))) + 1))
}
### source
mat <- mat[c("VenECs" , "EndoMA" , "CapECs" , "LECs" , "ArtECs") , 
           c("Proliferative Mac" , "FABP4+ Mac" , "SPP1+ Mac" , "SFTPC+ Mac" , "Monocytes" , "C9+ Mac" ,
             "infFib" , "apFib" , "myoFib" , "ecmFib" ,
             "Tm" , "yT" , "Proliferative T" , "NK" , "Anergy" , "Treg" , "Tn" , 
             "Plasma" , "FO-B")]
color.row <- c("VenECs" = "#F9EEBB" , "EndoMA" = "#C3E9F4" , "CapECs" = "#E5EAAC" , "LECs" = "#C9F2E6" , "ArtECs" = "#F2CCCC")
color.col <- rep(c("#EEE8E8" , "#d1e0e2" , "#d3d8e2" , "#d7e6d7") , c(6 , 4 , 7 , 2))
names(color.col) <- colnames(mat)
cell_groups <- c(rep(1:4 , c(6 , 4 , 7 , 2)))
coldf <- data.frame(group = colnames(mat))
rownames(coldf) <- colnames(mat)
rowdf <- data.frame(group = rownames(mat))
rownames(rowdf) <- rownames(mat)
ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                            border = FALSE, gp = gpar(fill = color.row, col = color.row)), 
                    show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                border = FALSE, gp = gpar(fill = color.col, col = color.col)), 
                        show_annotation_name = FALSE)
if (sum(abs(mat) > 0) == 1) {
  color.heatmap.use = c("white", color.heatmap.use)
}else {
  mat[mat == 0] <- NA
}
ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
              name = legend.name, top_annotation = ha2, 
              right_annotation = ha1, cluster_rows = F, 
              cluster_columns = F, row_names_side = "left", 
              row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
              column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
              column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
              row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), 
              row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                              fontface = "plain"), title_position = "leftcenter-rot", 
                                                              border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                                              grid_width = unit(2, "mm")) , 
              column_split = cell_groups , border = T)
ht1




