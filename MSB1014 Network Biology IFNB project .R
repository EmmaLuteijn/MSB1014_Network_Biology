# October 2024
# AD single-cell RNA-seq + network biology
# Project MSB1014 Network Biology

######## INSTALL PACKAGES ########
packages_to_install <- c("harmony", "Seurat", "ifnb.SeuratData", "tidyverse", "ggplot2","MAST", "biomaRt", "devtools", "CSCORE", "igraph", "RCy3", "clusterProfiler", "org.Hs.eg.db", "gridExtra")

for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

install_github("ChangSuBiostats/CS-CORE")

file_path <- "/Users/emmaluteijn/Documenten/R_studio/MSB1013/Project/Results/" # Set your own path
setwd(file_path) # Set the directory to this path

######## GET DATA ########
data("ifnb") # This will load the infb dataset from the ifnb.SeuratData package
str(ifnb) # Shows the structure of the data

######## QC and FILTERING ########
# Get the percentage of reads for mitochondrial genes. High mitochondrial gene expression, can indicate damaged cells.
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-') 

# Filter data on total RNA counts (UMI), detected genes (features), and mitochondrial gene expression.
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)

# Standard Seurat pipeline to prepare for dimensionality reduction and visualization
ifnb.filtered <- NormalizeData(ifnb.filtered) # Normalizes the data
ifnb.filtered <- FindVariableFeatures(ifnb.filtered) # Finds most variable genes across cells
ifnb.filtered <- ScaleData(ifnb.filtered) # Scales and center the data
ifnb.filtered <- RunPCA(ifnb.filtered) # Runs PCA
ElbowPlot(ifnb.filtered) # Helps to visualize each PC to variance

# Runs UMPA. Based on the elbow plot, it uses the first 20 principal components. 
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca') 

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim') + ggtitle("UMAP Before Harmony")# Shows the UMAP before Harmony

######## HARMONY ########
# Harmony integrates the data based on the condition (STIM or CTRL), to remove batch effects
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions # Check if the results are properly stored

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony") # Extracts the Harmony embeddings
ifnb.harmony.embed[1:10,1:10] # Views the first 10 cells and their first 10 Harmony dimensions

# Do UMAP and clustering using Harmony embeddings instead of PCA
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>% # Runs UMAP, but now using the Harmony embeddings
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% # Finds nearest neighbors of cells, used to find clusters
  FindClusters(resolution = 0.5) # Performs clustering, to later find the cell types

after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim') + ggtitle("UMAP After Harmony") # Shows UMAP after Harmony

# Side-by-side comparison of UMAP, to show that the data is correctly integrated
plt_UMAP_before_after <- before|after
ggsave(paste0(file_path,"plt_UMAP_before_after.png"), plot = plt_UMAP_before_after, width = 16, height = 8, dpi = 300)

######## Find cell types ########
# visualize data
clusters_plot <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition_plot <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

condition_plot|clusters_plot # show side-by-side comparision of UMAP visualized with the clusters or the condition of the group

# The cells already have annotations provided in the metadata
View(ifnb.harmony@meta.data)

Idents(ifnb.harmony) <- ifnb.harmony@meta.data$seurat_annotations # Setting Idents (cell identities) as Seurat annotations provided
plt_UMAP_conditions <- DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE) + ggtitle("UMAP")# Visualize UMAP with the cell types
ggsave(paste0(file_path, "plt_UMAP_conditions.png") , plot = plt_UMAP_conditions, width = 10, height = 8, dpi = 300)

######## Find Markers for B cells ########
ifnb.harmony$celltype.cnd <- paste0(ifnb.harmony$seurat_annotations,'_', ifnb.harmony$stim) # Paste condition and cell type together
Idents(ifnb.harmony) <- ifnb.harmony$celltype.cnd # Sets cell identities as the combination of condition and cell type

DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE) # Visualizes UMAP with the new idents

# find markers STIM vs CTRL B CELLS
B_cells <- subset(ifnb.harmony, seurat_annotations == "B") # Filters only the B cells

# Selects differentially expressed genes between stimulated and control B Cells. It uses the MAST test to perform the differential expression analysis. 
B_markers_MAST <- FindMarkers(ifnb.harmony, ident.1 = 'B_STIM', ident.2 = 'B_CTRL', test.use = 'MAST') 

######## Expression data for STIM B cells ########
B_cells_STIM <- subset(ifnb.harmony, celltype.cnd == "B_STIM") # Finds only the stimulated B cells
expr_data_STIM <- GetAssayData(B_cells_STIM, slot = "data") # Extract the normalized gene expression data
expr_data_STIM <- as.data.frame(expr_data_STIM) 
expr_data_STIM <- expr_data_STIM[rownames(expr_data_STIM) %in% rownames(B_markers_MAST),] # Selects only the differentially expressed genes in stimulated B cells

######## Expression data for CTRL B cells ########
B_cells_CTRL <- subset(ifnb.harmony, celltype.cnd == "B_CTRL") # Finds only the control B cells
expr_data_CTRL <- GetAssayData(B_cells_CTRL, slot = "data") # Extract the normalized gene expression data
expr_data_CTRL <- as.data.frame(expr_data_CTRL)
expr_data_CTRL <- expr_data_CTRL[rownames(expr_data_CTRL) %in% rownames(B_markers_MAST),] # Selects only the differentially expressed genes in stimulated B cells

# Make network B cells STIM ---------------------
# Do CSCORE analysis
CSCORE_result_B_STIM <- CSCORE(B_cells_STIM, genes = rownames(B_markers_MAST))

# Obtain CS-CORE co-expression estimates
CSCORE_coexp_B_STIM <- CSCORE_result_B_STIM$est 

# Obtain BH-adjusted p values
CSCORE_p_B_STIM <- CSCORE_result_B_STIM$p_value # Extracts the p-values
p_matrix_BH_B_STIM <- matrix(0, length(rownames(B_markers_MAST)), length(rownames(B_markers_MAST))) # Creates a matrix of zeros
p_matrix_BH_B_STIM[upper.tri(p_matrix_BH_B_STIM)] = p.adjust(CSCORE_p_B_STIM[upper.tri(CSCORE_p_B_STIM)], method = "BH") # Fills the upper triangle of the matrix of zeros with adjusted p-values
p_matrix_BH_B_STIM <- p_matrix_BH_B_STIM + t(p_matrix_BH_B_STIM) # This ensures that the matrix is symmetric

# Set co-expression entries with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp_B_STIM[p_matrix_BH_B_STIM > 0.05] <- 0

graph_B_STIM <- igraph::graph_from_adjacency_matrix(as.matrix(CSCORE_coexp_B_STIM), weighted = TRUE, mode = "upper", diag = FALSE)

avg_expression_B_STIM <- rowMeans(expr_data_STIM)
V(graph_B_STIM)$avg_expression_B_STIM <- avg_expression_B_STIM[V(graph_B_STIM)$name]
createNetworkFromIgraph(graph_B_STIM, title="Correlation_Network_B_Cells_STIM", collection="Correlation", nodeAttrList = list("avg_expression_B_STIM"))

# Analyze the network
analyzeNetwork()
degree_B_STIM <- getTableColumns("node", columns = c("name", "Degree", "BetweennessCentrality"), network = "Correlation_Network_B_Cells_STIM")

top10_degree_STIM <- degree_B_STIM %>%
  arrange(desc(Degree)) %>%
  slice_head(n = 10)

top10_betweenness_STIM <- degree_B_STIM %>%
  arrange(desc(BetweennessCentrality)) %>%
  slice_head(n = 10)

hubgenes_STIM <- intersect(top10_degree_STIM$name, top10_betweenness_STIM$name)

######## Analysis clusters B cells STIM ########
# gene info
ensemble <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
# Retrieve the clusters from cytoscape
Clusters_B_STIM <- getTableColumns("node", columns = c("name", "__glayCluster","avg_expression_B_STIM"), network = "Correlation_Network_B_Cells_STIM--clustered")

# Filter out the two clusters that contain the most genes with an above average gene expression
threshold_STIM <- quantile(Clusters_B_STIM$avg_expression_B_STIM, 0.75)

# Identify clusters with the highest number of genes exceeding the threshold <- maybe use this
top_clusters_STIM <- Clusters_B_STIM %>%
  filter(avg_expression_B_STIM > threshold_STIM) %>%  # Filter for high gene expression
  group_by(`__glayCluster`) %>%             # Group by cluster
  summarise(gene_count = n()) %>%          # Count the number of genes in each cluster
  arrange(desc(gene_count)) %>%             # Sort by gene count in descending order
  slice_head(n = 2)                         # Select the top 2 clusters


print(top_clusters_STIM)

######## Analysis FIRST cluster B cells STIM ########
Cluster_1_B_STIM <- Clusters_B_STIM %>%
  filter(Clusters_B_STIM$`__glayCluster` == top_clusters_STIM$`__glayCluster`[1])

print(Cluster_1_B_STIM)

# Get the gene info for the genes in B cells STIM
add_gene_info_STIM <- getBM(
 attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id", "go_id", "name_1006", "namespace_1003"),
 filters = "external_gene_name",
 values = Clusters_B_STIM$name,
 mart = ensemble)

# Add the gene information to the genes in cluster 1 of B Cells STIM
Cluster_1_gene_info_STIM <- merge(Cluster_1_B_STIM, add_gene_info_STIM, 
                             by.x = "name", by.y = "external_gene_name", 
                             all.x = TRUE)
  
print(Cluster_1_gene_info_STIM)

# GO analysis cluster 1 of B Cells STIM
go_1_STIM <- enrichGO(gene = Cluster_1_gene_info_STIM$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)

# Convert to data frame
go_1_STIM_df <- as.data.frame(go_1_STIM)

# Filter for significant GO terms
sig_go_1_STIM <- go_1_STIM_df %>%
  filter(p.adjust < 0.05)

# View or save significant GO terms cluster 1 of B Cells STIM
print(sig_go_1_STIM)

plt_bar_go_1_STIM <- barplot(go_1_STIM, showCategory = 10) + ggtitle("GO Community 4 - STIM")
plt_dot_go_1_STIM <- dotplot(go_1_STIM, showCategory = 10) + ggtitle("GO Community 4 - STIM")

######## Analysis SECOND cluster B cells STIM ########
Cluster_2_B_STIM <- Clusters_B_STIM %>%
  filter(Clusters_B_STIM$`__glayCluster` == top_clusters_STIM$`__glayCluster`[2])

print(Cluster_2_B_STIM)

# Add the gene information to the genes in cluster 1 of B Cells STIM
Cluster_2_gene_info_STIM <- merge(Cluster_2_B_STIM, add_gene_info_STIM, 
                                  by.x = "name", by.y = "external_gene_name", 
                                  all.x = TRUE)

print(Cluster_2_gene_info_STIM)

# GO analysis cluster 1 of B Cells STIM
go_2_STIM <- enrichGO(gene = Cluster_2_gene_info_STIM$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)

# Convert to data frame
go_2_STIM_df <- as.data.frame(go_2_STIM)

# Filter for significant GO terms
sig_go_2_STIM <- go_2_STIM_df %>%
  filter(p.adjust < 0.05)

# View or save significant GO terms cluster 1 of B Cells STIM
print(sig_go_2_STIM)

plt_bar_go_2_STIM <- barplot(go_2_STIM, showCategory = 10) + ggtitle("GO Community 3 - STIM")
plt_dot_go_2_STIM <- dotplot(go_2_STIM, showCategory = 10) + ggtitle("GO Community 3 - STIM")

plt_dot_go_1_STIM | plt_dot_go_2_STIM
# --------------------------------------------------------------------------------

######## Visual style of the network STIM ########
#Visual style of the overview graph
defaults <- list(
  NODE_SHAPE = "ellipse",  # Default node shape
  NODE_SIZE = 40,          # Default node size
  EDGE_TRANSPARENCY = 120,  # Default edge transparency
  NODE_COLOR = "darkgrey",  # Node label color
  NODE_BORDER_WIDTH = 2.0
)

node_color_map_STIM <- mapVisualProperty(
  visual.prop = "Node Fill Color",
  table.column = "__glayCluster",
  mapping.type = "d",
  table.column.values = c(top_clusters_STIM$`__glayCluster`[1], top_clusters_STIM$`__glayCluster`[2]),
  visual.prop.values = c("violet", "cyan")
)

node_size_map_STIM <- mapVisualProperty(
  visual.prop = "Node Size",
  table.column = "Degree",
  mapping.type = "continuous",
  table.column.values = c(min(degree_B_STIM$Degree), max(as.data.frame(degree_B_STIM$Degree))),
  visual.prop.values = c("40", "100")
)

hubgenes_map_border_STIM <- mapVisualProperty(
  visual.prop = "Node shape",
  table.column = "id",
  mapping.type = "d",
  table.column.values = c("ISG15","ISG20","CYBA"),
  visual.prop.values = c("diamond")
)

createVisualStyle("Single Cell co-expression STIM", defaults = defaults, mappings = list(node_color_map_STIM, node_size_map_STIM, hubgenes_map_border_STIM))

setVisualStyle("Single Cell co-expression STIM", network = "Correlation_Network_B_Cells_STIM")

# visual style of the clustered network
defaults <- list(
  NODE_SHAPE = "ellipse",  # Default node shape
  NODE_SIZE = 40,          # Default node size
  EDGE_TRANSPARENCY = 120,  # Default edge transparency
  NODE_COLOR = "darkgrey",  # Node label color
  NODE_BORDER_WIDTH = 2.0
)

node_color_map_STIM_cl <- mapVisualProperty(
  visual.prop = "Node Fill Color",
  table.column = "avg_expression_B_STIM",
  mapping.type = "c",
  table.column.values = c(min(as.data.frame(avg_expression_B_STIM)), max(as.data.frame(avg_expression_B_STIM))),
  visual.prop.values = c("white", "red")
)

createVisualStyle("Single Cell co-expression STIM communities", defaults = defaults, mappings = list(node_color_map_STIM_cl))

setVisualStyle("Single Cell co-expression STIM communities", network = "Correlation_Network_B_Cells_STIM--clustered")

######## Create Network B Cells CTRL ########
CSCORE_result_B_CTRL <- CSCORE(B_cells_CTRL, genes = rownames(B_markers_MAST))

# Obtain CS-CORE co-expression estimates
CSCORE_coexp_B_CTRL <- CSCORE_result_B_CTRL$est

# Obtain BH-adjusted p values
CSCORE_p_B_CTRL <- CSCORE_result_B_CTRL$p_value
p_matrix_BH_B_CTRL = matrix(0, length(rownames(B_markers_MAST)), length(rownames(B_markers_MAST)))
p_matrix_BH_B_CTRL[upper.tri(p_matrix_BH_B_CTRL)] = p.adjust(CSCORE_p_B_CTRL[upper.tri(CSCORE_p_B_CTRL)], method = "BH")
p_matrix_BH_B_CTRL <- p_matrix_BH_B_CTRL + t(p_matrix_BH_B_CTRL)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp_B_CTRL[p_matrix_BH_B_CTRL > 0.05] <- 0

graph_B_CTRL <- igraph::graph_from_adjacency_matrix(as.matrix(CSCORE_coexp_B_CTRL), weighted=TRUE, mode = "upper", diag = FALSE )

avg_expression_B_CTRL <- rowMeans(expr_data_CTRL)
V(graph_B_CTRL)$avg_expression_B_CTRL <- avg_expression_B_CTRL[V(graph_B_CTRL)$name]
createNetworkFromIgraph(graph_B_CTRL, title="Correlation_Network_B_Cells_CTRL", collection="Correlation", nodeAttrList = list("avg_expression_B_CTRL"))

# Analyze the network
analyzeNetwork()
degree_B_CTRL <- getTableColumns("node", columns = c("name", "Degree"), network = "Correlation_Network_B_Cells_CTRL")

degree_B_CTRL <- getTableColumns("node", columns = c("name", "Degree", "BetweennessCentrality"), network = "Correlation_Network_B_Cells_CTRL")

top10_degree_CTRL <- degree_B_CTRL %>%
  arrange(desc(Degree)) %>%
  slice_head(n = 10)

top10_betweenness_CTRL <- degree_B_CTRL %>%
  arrange(desc(BetweennessCentrality)) %>%
  slice_head(n = 10)

hubgenes_CTRL <- intersect(top10_degree_CTRL$name, top10_betweenness_CTRL$name)


# Visual style of the graph
node_color_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Fill Color",
  table.column = "avg_expression_B_CTRL",
  mapping.type = "continuous",
  table.column.values = c(min(as.data.frame(avg_expression_B_CTRL)), max(as.data.frame(avg_expression_B_CTRL))),
  visual.prop.values = c("white", "red")
)

node_size_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Size",
  table.column = "Degree",
  mapping.type = "continuous",
  table.column.values = c(min(degree_B_CTRL$Degree), max(as.data.frame(degree_B_CTRL$Degree))),
  visual.prop.values = c("40", "100")
)

node_label_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Label",
  table.column = "name",
  mapping.type = "p",
  visual.prop.values = c("black")
)

hubgenes_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Border Paint",
  table.column = "id",
  mapping.type = "d",
  table.column.values = c("IGJ","RPL13A","LMNA","OAS2"),
  visual.prop.values = c("black")
)

hubgenes_map_border_CTRL <- mapVisualProperty(
  visual.prop = "Node Border Width",
  table.column = "id",
  mapping.type = "d",
  table.column.values = c("IGJ","RPL13A","LMNA","OAS2"),
  visual.prop.values = c("10")
)

createVisualStyle("Single Cell Correlation CTRL", defaults = defaults, mappings = list(node_label_map_CTRL, node_color_map_CTRL, node_size_map_CTRL, hubgenes_map_CTRL))

setVisualStyle("Single Cell Correlation CTRL", network = "Correlation_Network_B_Cells_CTRL")


######## Analysis clusters B cells CTRL ########
# Retrieve the clusters from cytoscape
Clusters_B_CTRL <- getTableColumns("node", columns = c("name", "__glayCluster","avg_expression_B_CTRL"), network = "Correlation_Network_B_Cells_CTRL--clustered")

# Filter out the two clusters that contain the most genes with an above average gene expression
threshold_CTRL <- quantile(Clusters_B_CTRL$avg_expression_B_CTRL, 0.75)

# Identify clusters with the highest number of genes exceeding the threshold
top_clusters_CTRL <- Clusters_B_CTRL %>%
  group_by(`__glayCluster`) %>%  # Group by cluster first
  summarise(genes_above_threshold = sum(avg_expression_B_CTRL > threshold_CTRL)) %>%  # Count genes exceeding threshold in each cluster
  arrange(desc(genes_above_threshold)) %>%  # Sort clusters by the number of genes exceeding the threshold
  slice_head(n = 2)  # Select the top 2 clusters

print(top_clusters_CTRL)

######## Analysis FIRST cluster B cells CTRL ########
Cluster_1_B_CTRL <- Clusters_B_CTRL %>%
  filter(Clusters_B_CTRL$`__glayCluster` == top_clusters_CTRL$`__glayCluster`[1])

print(Cluster_1_B_CTRL)

# Get the gene info for the genes in B cells CTRL
add_gene_info_CTRL <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id", "go_id", "name_1006", "namespace_1003"),
  filters = "external_gene_name",
  values = Clusters_B_CTRL$name,
  mart = ensemble)

# Add the gene information to the genes in cluster 1 of B Cells CTRL
Cluster_1_gene_info_CTRL <- merge(Cluster_1_B_CTRL, add_gene_info_CTRL, 
                                  by.x = "name", by.y = "external_gene_name", 
                                  all.x = TRUE)

print(Cluster_1_gene_info_CTRL)

# GO analysis cluster 1 of B Cells CTRL
go_1_CTRL <- enrichGO(gene = Cluster_1_gene_info_CTRL$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)

# Convert to data frame
go_1_CTRL_df <- as.data.frame(go_1_CTRL)

# Filter for significant GO terms
sig_go_1_CTRL <- go_1_CTRL_df %>%
  filter(p.adjust < 0.05)

# View or save significant GO terms cluster 1 of B Cells CTRL
print(sig_go_1_CTRL)

plt_bar_go_1_CTRL <- barplot(go_1_CTRL, showCategory = 10) + ggtitle("GO Community 14 - CTRL")
plt_dot_go_1_CTRL <- dotplot(go_1_CTRL, showCategory = 10) + ggtitle("GO Community 14 - CTRL")

######## Analysis SECOND cluster B cells CTRL ########
Cluster_2_B_CTRL <- Clusters_B_CTRL %>%
  filter(Clusters_B_CTRL$`__glayCluster` == top_clusters_CTRL$`__glayCluster`[2])

print(Cluster_2_B_CTRL)

# Add the gene information to the genes in cluster 2 of B Cells CTRL
Cluster_2_gene_info_CTRL <- merge(Cluster_2_B_CTRL, add_gene_info_CTRL, 
                                  by.x = "name", by.y = "external_gene_name", 
                                  all.x = TRUE)

print(Cluster_2_gene_info_CTRL)

# GO analysis cluster 2 of B Cells CTRL
go_2_CTRL <- enrichGO(gene = Cluster_2_gene_info_CTRL$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)

# Convert to data frame
go_2_CTRL_df <- as.data.frame(go_2_CTRL)

# Filter for significant GO terms
sig_go_2_CTRL <- go_2_CTRL_df %>%
  filter(p.adjust < 0.05)

# View or save significant GO terms cluster 2 of B Cells CTRL
print(sig_go_2_CTRL)

plt_bar_go_2_CTRL <- barplot(go_2_CTRL, showCategory = 10) + ggtitle("GO Community 1 - CTRL")
plt_dot_go_2_CTRL <- dotplot(go_2_CTRL, showCategory = 10) + ggtitle("GO Community 1 - CTRL")

plt_dot_go_1_CTRL | plt_dot_go_2_CTRL

######## Visual style of the network ########
#Visual style of the overview graph
defaults <- list(
  NODE_SHAPE = "ellipse",  # Default node shape
  NODE_SIZE = 40,          # Default node size
  EDGE_TRANSPARENCY = 120,  # Default edge transparency
  NODE_COLOR = "darkgrey",  # Node label color
  NODE_BORDER_WIDTH = 2.0
)

node_color_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Fill Color",
  table.column = "__glayCluster",
  mapping.type = "d",
  table.column.values = c(top_clusters_CTRL$`__glayCluster`[1], top_clusters_CTRL$`__glayCluster`[2]),
  visual.prop.values = c("violet", "cyan")
)

node_size_map_CTRL <- mapVisualProperty(
  visual.prop = "Node Size",
  table.column = "Degree",
  mapping.type = "continuous",
  table.column.values = c(min(degree_B_CTRL$Degree), max(as.data.frame(degree_B_CTRL$Degree))),
  visual.prop.values = c("40", "100")
)

hubgenes_map_border_CTRL <- mapVisualProperty(
  visual.prop = "Node shape",
  table.column = "id",
  mapping.type = "d",
  table.column.values = c("ISG15","ISG20","CYBA"),
  visual.prop.values = c("diamond")
)

createVisualStyle("Single Cell co-expression CTRL", defaults = defaults, mappings = list(node_color_map_CTRL, node_size_map_CTRL, hubgenes_map_border_CTRL))

setVisualStyle("Single Cell co-expression CTRL", network = "Correlation_Network_B_Cells_CTRL")

# visual style of the clustered network
defaults <- list(
  NODE_SHAPE = "ellipse",  # Default node shape
  NODE_SIZE = 40,          # Default node size
  EDGE_TRANSPARENCY = 120,  # Default edge transparency
  NODE_COLOR = "darkgrey",  # Node label color
  NODE_BORDER_WIDTH = 2.0
)

node_color_map_CTRL_cl <- mapVisualProperty(
  visual.prop = "Node Fill Color",
  table.column = "avg_expression_B_CTRL",
  mapping.type = "c",
  table.column.values = c(min(as.data.frame(avg_expression_B_CTRL)), max(as.data.frame(avg_expression_B_CTRL))),
  visual.prop.values = c("white", "red")
)

createVisualStyle("Single Cell co-expression CTRL communities", defaults = defaults, mappings = list(node_color_map_CTRL_cl))

setVisualStyle("Single Cell co-expression CTRL communities", network = "Correlation_Network_B_Cells_CTRL--clustered")
# --------------------------------------------------------------------------------

# Visualization 
all_dot_plots <-(plt_dot_go_1_STIM | plt_dot_go_1_CTRL)/ 
(plt_dot_go_2_STIM | plt_dot_go_2_CTRL)
ggsave(paste0(file_path, "all_dot_plot.png"), plot = all_dot_plots, width = 20, height = 15, dpi = 300)

all_bar_plots <- (plt_bar_go_1_STIM | plt_bar_go_2_STIM)/ 
(plt_bar_go_1_CTRL  | plt_bar_go_2_CTRL)
ggsave(paste0(file_path,"all_bar_plot.png"), plot = all_bar_plots, width = 20, height = 15, dpi = 300)
