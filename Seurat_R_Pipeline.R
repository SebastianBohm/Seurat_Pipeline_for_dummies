

#################################################
#################################################
#### SEBASTIAN's EASY TO USE SEURAT PIPELINE #### :)
#################################################
#################################################

# If questions/problems ask ChatGPT or me

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

# Define output and working directories
output_dir <- "XXX/"
work_path <- "XXX/"

# Load Seurat object
seurat_obj <- readRDS("XXX/seurat_obj_combined_v2.rds")

#######################################################
#
### STANDARD ANALYSIS (do not change)
# [You have to run this, before running something else]
#######################################################

# Visualize t-SNE and UMAP
DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + DimPlot(seurat_obj, reduction = "umap", label = TRUE)

umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
ggsave(paste0(output_dir, "/AAA_seurat_UMAP_cluster.pdf"), plot = umap_plot, device = "pdf", width = 7, height = 5)

# UMAP grouped by RT_group
umap_time_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "RT_group", label = FALSE)
ggsave(paste0(output_dir, "/AAA_seurat_UMAP_time.pdf"), plot = umap_time_plot, device = "pdf", width = 7, height = 5)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.9)
umap_plot_0.9 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Resolution 0.9")
ggsave(paste0(output_dir, "/seurat_UMAP_res_0.9.jpg"), plot = umap_plot_0.9, device = "jpg", width = 7, height = 5)
umap_plot_0.9
#markers_0.9 <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top_markers_0.9 <- markers_0.9 %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
#write.csv(top_markers_0.9, paste0(output_dir, "/seurat_top_25_markers_res_0.9.csv"), row.names = FALSE)




#########################################################################
##
## RUN HEATMAPS of Genes of interest (change "Genes")
## [Look up gene names in "annotation_file"]
##
annotation_file <- "XXX/df_gene.csv"
gene_annotation <- read.csv(annotation_file)
########################################################################

# Define the list of genes
#############################
############################# ONLY CHANGE THIS
#############################
genes <- c("BMP5", "ALX1", "VEGFA", "TBX5","TBX5[nr]", "CDH5","LOC104376251[nr]|CDH5[hs]", "BMPER[hs]", "LOC106489042[nr]|BMPER[hs]", "DND1", "LOC112105876[nr]|AXDND1[hs]", "DAZL", "CD34", "KDR")

# Assuming your annotation file contains columns "gene_short_name" and "gene_id" (with AMEX60DD IDs)
# Filter the annotation table to only keep relevant genes
genes_filtered <- gene_annotation %>% filter(gene_short_name %in% genes)

# Match the gene IDs (AMEX60DD) with the rownames in your Seurat object
genes_found <- intersect(rownames(seurat_obj@assays$RNA), genes_filtered$gene_id)

# Loop through each matched gene and create the plots
for (gene_id in genes_found) {
  # Get the gene short name for the plot title
  gene_short_name <- genes_filtered$gene_short_name[genes_filtered$gene_id == gene_id]
  
  # Check if the gene has any non-zero expression
  gene_data <- FetchData(seurat_obj, vars = gene_id)
  
  if (all(gene_data == 0)) {
    message(paste("Skipping", gene_short_name, "due to no expression."))
    next  # Skip this gene and continue with the next one
  }
  
  # UMAP Feature Plot (heatmap-style plot)
  p1 <- FeaturePlot(seurat_obj, features = gene_id, reduction = "umap", cols = c("lightgrey", "blue")) +
    ggtitle(paste0("UMAP Feature Plot: ", gene_short_name))
  ggsave(filename = paste0(output_dir, gene_short_name, "_heat_umap.jpg"), plot = p1, device = "jpg")
  
  # Violin plot for the gene expression
  p2 <- VlnPlot(seurat_obj, features = gene_id, pt.size = 0) +
    ggtitle(paste0("Violin Plot: ", gene_short_name))
  ggsave(filename = paste0(output_dir, gene_short_name, "_violin.jpg"), plot = p2, device = "jpg")
}


##############################################################################################
##
## SUB-CLUSTERING
## [change "idents" to you clusters of interest]
##############################################################################################

# Subset the Seurat object to include clusters 1, 2, 10, 13, 17, 18, 19, 20, 21, 24, 27
subcluster_custom <- subset(seurat_obj, idents = c(1, 2, 10, 13, 17, 18, 19, 20, 21, 24, 27))

# Re-run clustering for this subset
subcluster_custom <- FindVariableFeatures(subcluster_custom)
subcluster_custom <- ScaleData(subcluster_custom)
subcluster_custom <- RunPCA(subcluster_custom)
subcluster_custom <- RunUMAP(subcluster_custom, dims = 1:10)
subcluster_custom <- FindNeighbors(subcluster_custom, dims = 1:10)
subcluster_custom <- FindClusters(subcluster_custom, resolution = 0.5)

# Save the initial UMAP without group.by
initial_umap_custom <- DimPlot(subcluster_custom, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Initial UMAP: Clusters 1, 2, 10, 13, 17, 18, 19, 20, 21, 24, 27")
ggsave("XXXX/initial_umap_subcluster_custom.jpg", plot = initial_umap_custom, device = "jpg", width = 7, height = 5)

# Generate the UMAP plot grouped by RT_group
umap_plot_custom <- DimPlot(subcluster_custom, reduction = "umap", group.by = "RT_group", repel = TRUE) + 
  ggtitle("Subcluster: Clusters 1, 2, 10, 13, 17, 18, 19, 20, 21, 24, 27")
ggsave("XXX/umap_subcluster_custom.jpg", plot = umap_plot_custom, device = "jpg", width = 7, height = 5)

# Save the normal UMAP with cluster labels
umap_normal_custom <- DimPlot(subcluster_custom, reduction = "umap", label = TRUE) + 
  ggtitle("Resolution: Subcluster 1, 2, 10, 13, 17, 18, 19, 20, 21, 24, 27")
ggsave("XXX/normal_umap_subcluster_custom.jpg", plot = umap_normal_custom, device = "jpg", width = 7, height = 5)

# Find top 25 marker genes for each subcluster
markers_custom <- FindAllMarkers(subcluster_custom, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers_custom <- markers_custom %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)

# Join with gene names
top_markers_custom_with_names <- top_markers_custom %>%
  left_join(df_gene, by = c("gene" = "gene_id")) %>%
  select(-gene_type, -chr)

write.csv(top_markers_custom_with_names, paste0(output_dir, "/top_25_marker_genes_subcluster_custom.csv"), row.names = FALSE)

####################################
###################################
####################################



