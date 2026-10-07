library(nichenetr)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(paletteer)
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS('lr_network.rds')
weighted_networks = readRDS('weighted_networks.rds')

data <- readRDS("supercell.rds")
### ecmFib
sender_celltypes <- unique(data@meta.data$metacell_type)[-c(6,20,21,24)] ### "undifined"
nichenet_output_ecmFib <- nichenet_seuratobj_aggregate(
  seurat_obj = data , 
  receiver = c("ecmFib") , 
  condition_colname = "disease" , condition_oi = "Fibrosis" , condition_reference = "Normal" , 
  sender = sender_celltypes , 
  ligand_target_matrix = ligand_target_matrix , lr_network = lr_network , weighted_networks = weighted_networks)

nichenet_output_ecmFib$ligand_activity_target_heatmap

saveRDS(nichenet_output_ecmFib , "ecmFib_nichenet.rds")

