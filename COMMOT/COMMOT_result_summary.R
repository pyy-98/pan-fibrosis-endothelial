library(Seurat)
library(CellChat)
data <- readRDS("./Niche_expression.rds")

CellChatDB <- CellChatDB.human
cellchat <- CellChatDB[["interaction"]]

## example
result <- list()
for (i in c("V19S23-092-A1" , "V19S23-092-B1" , "V19S23-092-C1" , "V19S23-092-D1")) {
  COMMOT <- read.table(paste("./COMMOT/" , i , "_receiver_matrix.txt" , sep = "") , sep = "\t" , header = T , row.names = 1)  ## sender
  subdata <- subset(data , subset = sample_id == i)
  subdata <- subset(subdata , idents = "niche_0")
  COMMOT <- COMMOT[unlist(lapply(strsplit(colnames(subdata) , "_") , "[" , 2)),]
  colnames(COMMOT) <- gsub("r\\." ,"", colnames(COMMOT)) ## s\\.
  col_pathway_name <- unlist(lapply(strsplit(colnames(COMMOT) , "\\.") , "[" , 1))
  interact <- data.frame(pathway_name = unique(cellchat$pathway_name))
  interact$score <- rep(0 , nrow(interact))
  for (j in interact$pathway_name) {
    subcellchat <- cellchat[which(cellchat$pathway_name == j),]
    sub_COMMOT <- COMMOT[,col_pathway_name %in% subcellchat$ligand]
    interact$score[which(interact$pathway_name == j)] <- sum(sub_COMMOT)/nrow(COMMOT)
  }
  result[[i]] <- interact
}

final <- data.frame(Normal = result[[1]]$score , 
                    Mild = result[[2]]$score , 
                    Moderate = result[[3]]$score , 
                    Severe = result[[4]]$score , row.names = result[[1]]$pathway_name)



final <- final[order(rowMeans(final) , decreasing = T)[1:20],]
library(pheatmap)
library(RColorBrewer) 
annotation_col <- data.frame(row.names = colnames(final) , IPF_level = colnames(final))
ann_colors <- list(
  IPF_level = c(Normal = "#f3e0d6" , Mild = "#e7af97" , Moderate = "#c57466" , Severe = "#9c363b")
)

p <- pheatmap(final , 
              scale = "column" , 
              color = colorRampPalette(c("#6799A5" , "white" , "#C1687B"))(50) ,
              annotation_col = annotation_col ,  
              annotation_colors = ann_colors , 
              fontsize_col = 10 , 
              cluster_cols = F , 
              cluster_rows = F , 
              show_rownames = T , 
              show_colnames = F , 
              fontsize = 10 , 
              cellwidth = 10 , 
              cellheight = 10 , 
              border = F , 
              main = "Receiver"
)
p

