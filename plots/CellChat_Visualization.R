library(CellChat)
library(Seurat)
library(ggplot2)
cellchat_Fibrosis <- readRDS("Fibrosis_cellchat.rds")
cellchat_Normal <- readRDS("Normal_cellchat.rds")

object.list <- list(Normal = cellchat_Normal , Fibrosis = cellchat_Fibrosis)
cellchat <- mergeCellChat(object.list, add.names = names(object.list) , cell.prefix = T)

### 互作强度热图
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

#### 不同细胞类型之间相互作用或交互强度的差异
options(repr.plot.height = 6, repr.plot.width = 12)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)



### 展示两组间上下调的配受体对
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("COLLAGEN"))
netVisual_bubble(cellchat , 
                 sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.01 , pairLR.use = pairLR.use)


pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL"))
netVisual_bubble(cellchat , 
                 sources.use = c("C9+ Mac" , "Proliferative Mac" , "SPP1+ Mac" , "SFTPC+ Mac" , "FABP4+ Mac" , "Monocytes") , targets.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , 
                 comparison = c(1,2) , thresh = 0.05 , pairLR.use = pairLR.use)

netVisual_bubble(cellchat , 
                 sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.001)


pairLR.use <- extractEnrichedLR(cellchat, signaling = c("FN1"))
netVisual_bubble(cellchat , 
                 sources.use = c("C9+ Mac" , "Proliferative Mac" , "SPP1+ Mac" , "SFTPC+ Mac" , "FABP4+ Mac" , "Monocytes") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.05 , pairLR.use = pairLR.use)


netVisual_bubble(cellchat , 
                 sources.use = c("C9+ Mac" , "Proliferative Mac" , "SPP1+ Mac" , "SFTPC+ Mac" , "FABP4+ Mac" , "Monocytes") , targets.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , 
                 comparison = c(1,2) , thresh = 0.05)

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MHC-II"))
netVisual_bubble(cellchat , 
                 sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.05 , pairLR.use = pairLR.use)

netVisual_heatmap(cellchat , signaling = c("MHC-II") , measure = "weight" , 
                  sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , 
                  targets.use = c("myoFib") , title.name = "MHC-II signaling differential interaction strength")

plotGeneExpression(cellchat, features = c("CCL2" , "CCL5" , "CCL7" , "CCL8" , "CCL13", "CCL17" , "CCL18" , 
                                          "ACKR1" , "ACKR2") , 
                   color.use = color <- c("FO-B" = "#a9bfb7" , "Plasma" = "#65848a" , "infFib" = "#b97f52" , 
                                          "myoFib" = "#738a98" , "apFib" = "#8f493b" , "ecmFib" = "#706e92" ,  
                                          "Proliferative Mac" = "#8c8aa6" , "SPP1+ Mac" = "#4e7698" , "Monocytes" = "#7fa781" , 
                                          "C9+ Mac" = "#656565" , "SFTPC+ Mac" = "#7293ac" , "FABP4+ Mac" = "#5d8d61" , 
                                          "Treg" = "#bb6c70" , "NK" = "#dbca93" , "yT" = "#a1c3c1" , "Tn" = "#a9cad9" , "Tm" = "#5c6841" , 
                                          "Proliferative T" = "#966a93" , "Anergy" = "#c1875b" , 
                                          "VenECs" = "#F9EEBB" , "EndoMA" = "#C3E9F4" , "CapECs" = "#E5EAAC" , "LECs" = "#C9F2E6" , "ArtECs" = "#F2CCCC")
)


pairLR.use <- extractEnrichedLR(cellchat, signaling = c("FN1"))
netVisual_bubble(cellchat , 
                 sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.05 , pairLR.use = pairLR.use)

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("LAMININ"))
pairLR.use <- pairLR.use[c(grep("LAMB1",pairLR.use$interaction_name),grep("LAMB2",pairLR.use$interaction_name)),,drop = F]
netVisual_bubble(cellchat , 
                 sources.use = c("VenECs" , "ArtECs" , "CapECs" , "LECs" , "EndoMA") , targets.use = c("infFib" , "myoFib" , "apFib" , "ecmFib") , 
                 comparison = c(1,2) , thresh = 0.05 , pairLR.use = pairLR.use)



