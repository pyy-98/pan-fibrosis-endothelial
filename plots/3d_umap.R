library(scatterplot3d)
library(tidyverse)
library(openxlsx)
library(factoextra)
library(Seurat)
### data is seurat.rds of scRNAseq expression matrice of interest cells
data <- readRDS("...")
data <- RunUMAP(data, dims = 1:10,                
                reduction = "pca",                
                n.neighbors=10,                 
                min.dist=0.5,                
                n.components = 5)
UMAP_3d = data@reductions$umap@cell.embeddings %>%  as.data.frame() %>%  cbind(celltype=data@meta.data$metacell_type)
library(dplyr)
library(scatterplot3d)
color.bin <- c("FO-B" = "#a9bfb7" , "Plasma" = "#65848a" , "infFib" = "#b97f52" , 
               "myoFib" = "#738a98" , "apFib" = "#8f493b" , "ecmFib" = "#706e92" ,  
               "Proliferative Mac" = "#8c8aa6" , "SPP1+ Mac" = "#4e7698" , "Monocytes" = "#7fa781" , 
               "C9+ Mac" = "#656565" , "SFTPC+ Mac" = "#7293ac" , "FABP4+ Mac" = "#5d8d61" , 
               "Treg" = "#bb6c70" , "NK" = "#dbca93" , "yT" = "#a1c3c1" , "Tn" = "#a9cad9" , "Tm" = "#5c6841" , 
               "Proliferative T" = "#966a93" , "Anergy" = "#c1875b" , "VenECs" = "#F9EEBB" , "EndoMA" = "#C3E9F4" , 
               "CapECs" = "#E5EAAC" , "LECs" = "#C9F2E6" , "ArtECs" = "#F2CCCC")
color.bin <- data.frame(color = color.bin , celltype = names(color.bin))
mycol <- color.bin$color[match(UMAP_3d[,6] , color.bin$celltype)]

scatterplot3d(UMAP_3d , 
              color = mycol , pch = 16, cex.symbols = 0.1,              
              scale.y = 0.7,angle = 30,            
              xlab = "UMAP_2", ylab = "UMAP_1", zlab = "UMAP_3",              
              main="3D Scatter Plot of UMAP",              
              col.axis = "#444444", col.grid = "#CCCCCC")
legend("bottom", legend = levels(as.factor(UMAP_3d$celltype)),       
       col =  color.bin$color[match(levels(as.factor(UMAP_3d$celltype)) , color.bin$celltype)],  pch = 16,       
       inset = -1,xpd = TRUE, horiz = FALSE,ncol = 3)




