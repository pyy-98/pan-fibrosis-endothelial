library(Seurat)
library(SuperCell)
##  seurat_list is a list of scRNAsea expression matrices of different cell types
seurat_list <- readRDS("...")

data <- merge(x = seurat_list[[1]] , y = seurat_list[-1])
exp_mat <- GetAssayData(data , slot = "data")
data <- FindVariableFeatures(data , nfeatures = 2000)
hvg <- VariableFeatures(data)
SC <- SCimplify(exp_mat , k.knn = 5 , gamma = 10 , genes.use = hvg)
SC.GE <- supercell_GE(exp_mat , SC$membership)
SC$metacell_type <- supercell_assign(clusters = data@meta.data$Cell_type , 
                                 supercell_membership = SC$membership , 
                                 method = "jaccard")
purity <- supercell_purity(clusters = data@meta.data$Cell_type , supercell_membership = SC$membership , 
                           method = "entropy")

supercell <- supercell_2_Seurat(SC.GE = as.matrix(SC.GE) , SC = SC , 
                                         fields = c("metacell_type" , "purity"))

supercell <- RunUMAP(supercell , dims = 1:10)
supercell <- FindClusters(supercell , graph.name = "RNA_nn")
Idents(supercell) <- supercell@meta.data$metacell_type
color <- c("FO-B" = "#a9bfb7" , "Plasma" = "#65848a" , "infFib" = "#b97f52" , 
           "myoFib" = "#738a98" , "apFib" = "#8f493b" , "ecmFib" = "#706e92" ,  
           "Proliferative Mac" = "#8c8aa6" , "SPP1+ Mac" = "#4e7698" , "Monocytes" = "#7fa781" , 
           "C9+ Mac" = "#656565" , "SFTPC+ Mac" = "#7293ac" , "FABP4+ Mac" = "#5d8d61" , 
           "Treg" = "#bb6c70" , "NK" = "#dbca93" , "yT" = "#a1c3c1" , "Tn" = "#a9cad9" , "Tm" = "#5c6841" , 
           "Proliferative T" = "#966a93" , "Anergy" = "#c1875b" , "EndoMA" = "#F2C8E1" , "LECs" = "#C8E6C9" ,
           "ArtECs" = "#FCE4EC" , "VenECs" = "#FFE4B5" , "CapECs" = "#F7B7A3")
DimPlot(supercell , label = T , cols = color , raster = F)



