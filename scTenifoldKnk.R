library(scTenifoldKnk)
library(Seurat)
library(dplyr)
library(forcats)
library(ggplot2)
### data is scRNAseq expression matrices of interest cell type
data <- readRDS("...")
# data <- subset(data , idents = "VenECs")
# data <- subset(data , downsample = 1000)

countMatrix <- GetAssayData(data , slot = "counts")
# col <- dim(countMatrix)[2]
# Ei <- apply(countMatrix , 1, function(x){sum(as.numeric(x)/col)})
# countMatrix <- countMatrix[order(Ei , decreasing = T)[1:10000] , ]

result <- scTenifoldKnk(countMatrix = countMatrix , 
                        gKO = "Npm1" , ## gene of interest
                        nCores = 24)

top_genes <- result$diffRegulation[which(result$diffRegulation$p.adj < 0.05),]
top_genes <- top_genes[-grep("mt-" , top_genes$gene),]

top_genes$logFC <- log2(top_genes$FC)
top_genes$neg_log10_pval <- -log10(top_genes$p.adj)
top_genes$gene <- fct_reorder(top_genes$gene , top_genes$logFC)
top_genes[sapply(top_genes , is.infinite)] <- 260

p<-ggplot(top_genes, aes(x = logFC, y = gene)) +
  geom_bar(
    aes(fill = neg_log10_pval),
    stat = "identity",
    width = 0.7,
    alpha = 0.8
  ) +
  geom_text(
    aes(label = sprintf("%.2f", logFC)),
    hjust = -0.2,
    size = 3.5,
    fontface="bold"
  ) +
  
  #设置颜色
  scale_fill_gradient(
    low="lightblue",
    high="darkblue",
    name="-log10(P)",
    guide = guide_colorbar(
      title.position="top",
      title.hjust = 0.5
    )) +
  scale_y_discrete(limits = rev) +
  # 设置坐标轴和标题
  labs(
    x="Log2(Perturbation Magnitude)",
    y="Gene",
    title="Top Perturbed Genes" , 
    subtitle = "KO:Npm1"
  ) +
  # 自定义主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color="gray40"),
    axis.title = element_text(size = 14, face="bold"),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(face="italic"), # Italicize gene names
    legend.position="top",
    legend.title = element_text(size = 12, face="bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color="gray90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 40, 20, 20) # Adjust margins for labels
  ) + expand_limits(x = max(top_genes$logFC) * 1.2)
p
