library(Seurat)
library(ggplot2)
library(progeny)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(circlize)
data <- readRDS("./Endothelial.rds")
data <- progeny(data , scale = F , organism = "Human" , top = 500 , perm = 1 , return_assay = T)
data <- ScaleData(data , assay = "progeny")

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(data, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

CellsClusters <- data.frame(Cell = names(Idents(data)), 
                            CellType = as.character(data@meta.data$Cell_type),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df , CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
write.table(summarized_progeny_scores_df , file = "PROGENy.txt" , sep = "\t" , quote = F)



