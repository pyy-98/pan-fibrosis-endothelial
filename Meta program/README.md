# Source code for Gavish et al.

## Generating_MPs
- Generate_Meta_Programs.R generates MPs from NMFs programs that were calculated for each sample using different ranks. The NMF programs were calculated per sample using the ‘NMF’ R package:  
NMFs_per_sample = nmf(x = Expression_matrix, rank = 4:9, method="snmf/r", nrun = 10)  
NMF programs are listed in Genes_nmf_w_basis, where each entry contains NMF gene-scores of a single sample.  
Using the function robust_nmf_programs.R (described below) and a custom written clustering method, MPs are derived from the NMF programs.
- Genes_nmf_w_basis_example.RDS is an example for how the NMF output is arranged. Note that each element name in the list ends with '_rank4_9_nruns10.RDS', and each matrix column ends with an extension that represents the NMF rank and program index. 
- robust_nmf_programs.R filters outs redundant and non-robust NMF programs for each sample. The remaining NMF programs are those that appeared in more than one NMF rank within the sample and have similarity to a NMF program in a different sample. Please see also https://github.com/gabrielakinker/CCLE_heterogeneity for more details on how to define robust NMF programs.

## MPs_distribution
- MP_distribution.R depicts the meta-program distribution of a certain cell type within individual samples and across the whole study. It requires a list of scRNAseq expression matrices of different samples in a study (of a certain cell type).
- MP_list.RDS is a list of meta-programs for each cell type required for inferring the distributions in MP_distribution.R
- heatCols.RDS is the colormap used for plotting heatmaps in MP_distribution.R
- My_study.RDS can be used for an example
