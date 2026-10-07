library(Seurat)
library(CellChat)
supercell <- readRDS("supercell.rds")
Normal_data <- GetAssayData(supercell , assay = "RNA" , slot = "data")
Normal_meta <- data.frame(group = supercell@meta.data$metacell_type , row.names = colnames(Normal_data))

cellchat <- createCellChat(object = Normal_data  , meta = Normal_meta , group.by = "group")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat , slot.name = "netP")
saveRDS(cellchat , "cellchat.rds")





