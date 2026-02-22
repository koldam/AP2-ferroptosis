library(monocle3)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)

# Step 1: Load the data
expression_matrix = read.table("___.txt", sep = "\t", dec = ".", header = T)
rownames(expression_matrix) = expression_matrix$X
expression_matrix$X = NULL
cell_metadata = read.table("___.txt", sep = "\t", dec = ".", header = T)
rownames(cell_metadata) = cell_metadata$X
cell_metadata$X = NULL
gene_annotation = read.table("___.txt", sep = "\t", dec = ".", header = T)
rownames(gene_annotation) = gene_annotation$X
gene_annotation$X = NULL

library(Matrix)
matrix = as(exprs, "dgCMatrix")

cds = new_cell_data_set(matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 2: Normalize and pre-process the data
cds = preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

## Step 3: Reduce the dimensions using UMAP
cds = reduce_dimension(cds)
plot_cells(cds)

## Step 4: Cluster the cells
cds = cluster_cells(cds)

## Step 5: Analyze
pr_graph_test_res = graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df = find_gene_modules(cds[pr_deg_ids,], k=4)

agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c(" ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="row", clustering_method="ward.D2",
                   fontsize=8)
