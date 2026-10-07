library(Seurat)
library(hdf5r)
library(ggplot2)
library(spacexr)
library(readxl)
dir <- "./E-MTAB-14121/"
seurat_list <- readRDS("E-MTAB-14121.rds")
clinical <- read_xlsx("Clinical.xlsx")
clinical <- clinical[which(clinical$Source == "E-MTAB-14121"),]
clinical <- clinical[which(clinical$Disease == "Normal"),] # Fibrosis
seurat_list <- seurat_list[clinical$Sample]
### ref
ref <- readRDS("Normal_ref.rds") # Fibrosis_ref.rds
counts <- ref[["RNA"]]@counts
cluster <- Idents(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

for(i in 1:length(seurat_list)){
  ### obs
  data <- seurat_list[[i]]
  counts <- data[["Spatial"]]@counts
  rownames(counts) <- toupper(rownames(counts))
  coords <- GetTissueCoordinates(data)
  coords <- coords[,c(1,2)]
  colnames(coords) <- c("x" , "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  RCTD <- create.RCTD(query, reference, max_cores = 2)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  saveRDS(RCTD , file = paste(dir , "RCTD_result/" ,  names(seurat_list)[i] , "_RCTD_result.rds" , sep = ""))
  print(i)
}