# Run Harmony integration
library(gridExtra)
library(conflicted)
library(dplyr)
library(harmony)
library(reshape2)
sc_569 <- SCTransform(sc_569)
sc_569 <- RunPCA(sc_569, npcs = 30, verbose = F)
sc_569 <- FindNeighbors(sc_569, dims = 1:30, reduction = "integrated.dr")
sc_569 <- FindClusters(sc_569, resolution = 2)

harmony_embeddings <- HarmonyMatrix(
  data_mat = Embeddings(sc_569, 'pca'),
  meta_data = sc_569@meta.data,
  vars_use = "group",
  do_pca = FALSE
)
# Add harmony embeddings to Seurat sc_569ect
sc_569[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = "harmony_",
  assay = "SCT"
)

# Run UMAP using harmony reduction
sc_569 <- RunUMAP(sc_569, reduction = "pca", dims = 1:30, 
                  reduction.name = "umap")
# Find neighbors using harmony reduction
sc_569 <- FindNeighbors(sc_569, reduction = "pca", dims = 1:30)
saveRDS(seurat_obj, file = "GSE223128_sc569_harmonized.rds")
sc_569 <- FindClusters(sc_569, resolution = 0.2, algorithm = 4, verbose = FALSE)

# Create UMAP plot for this cluster
print(DimPlot(
  sc_569,
  reduction = "umap",
  group.by = "group",
  label = FALSE
) +
  ggtitle("UMAP by Replicate") +
  theme_minimal())

# Create main plot directory if it doesn't exist
main_plot_dir <- "/home/saini_lab/Documents/mitoDynamics/plots/GSE223128/Cluster visualizations"
if (!dir.exists(main_plot_dir)) {
  dir.create(main_plot_dir, recursive = TRUE)
}

resolutions <- seq(0.1, 0.9, by = 0.1)
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")
# Acquire the G1 phase genes
g1_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G1") %>%
  pull("gene_name")
# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

sc_569 <- CellCycleScoring(sc_569,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)

# Loop through different resolutions
for (res in resolutions) {
  res_folder <- paste0("res_", res)
  
  # Create directory for this resolution if it doesn't exist
  res_path <- file.path(main_plot_dir, res_folder)
  if (!dir.exists(res_path)) {
    dir.create(res_path, recursive = TRUE)
  }
  
  # Find clusters at current resolution
  cluster_col <- paste0("harmony_clusters_res", res)
  
  sc_569 <- FindClusters(sc_569, resolution = res, 
                         cluster.name = cluster_col)
  sc_569 <- CellCycleScoring(sc_569, g2m.features = g2m_genes,
                             s.features = s_genes)
  
  # Create UMAP plot with clusters
  pdf(file.path(res_path, paste0("sc_569_cluster_umap_res", res, ".pdf")), width = 10, height = 8)
  print(DimPlot(sc_569, reduction = "umap.harmony", 
                group.by = cluster_col, 
                label = TRUE, raster = FALSE) +
          ggtitle(paste0(sample_name, " - Clusters (Resolution ", res, ")")))
  dev.off()
  
  # Create feature plots for QC metrics
  pdf(file.path(res_path, paste0("sc_569_qc_metrics_res", res, ".pdf")), width = 10, height = 8)
  
  # Plot nFeature (genes detected)
  p1 <- FeaturePlot(sc_569, features = "nFeature_SCT", reduction = "umap.harmony") +
    ggtitle("Gene count")
  
  # Plot nCount (UMI count)
  p2 <- FeaturePlot(sc_569, features = "nCount_SCT", reduction = "umap.harmony") +
    ggtitle("UMI count")
  
  # Plot percent mitochondrial
  p3 <- FeaturePlot(sc_569, features = "mitoRatio", reduction = "umap.harmony") +
    ggtitle("% Mitochondrial")
  
  # Plot cell cycle scores and phase
  p4 <- DimPlot(sc_569, reduction = "umap.harmony", group.by = "Phase") +
    ggtitle("Cell Cycle Phase")
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  
  # Fix for proportion plots
  # Make sure sample column exists
  if (!"sample" %in% colnames(sc_569@meta.data)) {
    sc_569@meta.data$sample <- sample_name
  }
  
  metadata <- sc_569@meta.data
  
  # Create proportion plots PDF
  pdf(file.path(res_path, paste0("sc_569_proportion_plots_res", res, ".pdf")), width = 10, height = 8)
  
  # 1. Cluster distribution plot
  cluster_counts <- metadata %>%
    group_by(sample, cluster = .data[[cluster_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p_clusters <- ggplot(cluster_counts, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    theme_minimal() +
    labs(title = paste0(sample_name, " (res ", res, "): Cluster Distribution"),
         x = "Sample", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_clusters)
  
  # 2. Cell cycle distribution plot
  phase_counts <- metadata %>%
    group_by(sample, Phase) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p_phase <- ggplot(phase_counts, aes(x = sample, y = proportion, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    theme_minimal() +
    labs(title = paste0(sample_name, ": Cell Cycle Distribution"),
         x = "Sample", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_phase)
  
  # 3. Cell cycle by cluster plot (avoiding pivot_wider)
  phase_by_cluster <- metadata %>%
    group_by(cluster = .data[[cluster_col]], Phase) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p_phase_by_cluster <- ggplot(phase_by_cluster, aes(x = cluster, y = proportion, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    theme_minimal() +
    labs(title = paste0(sample_name, " (res ", res, "): Cell Cycle by Cluster"),
         x = "Cluster", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_phase_by_cluster)
  
  # Combined plot
  combined_plot <- grid.arrange(p_clusters, p_phase, ncol = 2,
                                top = textGrob(paste0(sample_name, " - Resolution ", res),
                                               gp = gpar(fontsize = 15, fontface = "bold")))
  print(combined_plot)
  
  dev.off()
  
  # Write number of clusters
  n_clusters <- length(unique(metadata[[cluster_col]]))
  writeLines(
    paste("Number of clusters at resolution", res, ":", n_clusters),
    file.path(res_path, paste0("sc_569_cluster_stats_res", res, ".txt"))
  )
}
# Generate additional sample-specific QC plots
pdf(file.path(main_plot_dir, "sc_569_qc_metrics_overview.pdf"), width = 12, height = 10)

# Create violin plots for QC metrics
VlnPlot(sc_569, features = c("nFeature_SCT", "nCount_SCT", "mitoRatio"), 
        ncol = 2, pt.size = 0) &
  theme(plot.title = element_text(size = 12, face = "bold"))

# Scatter plots for QC metrics correlations
p1 <- FeatureScatter(sc_569, feature1 = "nCount_SCT", feature2 = "nFeature_SCT") +
  ggtitle("Gene count vs UMI count")

p2 <- FeatureScatter(sc_569, feature1 = "nCount_SCT", feature2 = "mitoRatio") +
  ggtitle("UMI count vs % Mitochondrial")

p3 <- FeatureScatter(sc_569, feature1 = "nFeature_SCT", feature2 = "mitoRatio") +
  ggtitle("Gene count vs % Mitochondrial")

print(p1 + p2)
print(p3 + p4)

# Cell cycle score correlation
p5 <- FeatureScatter(sc_569, feature1 = "S.Score", feature2 = "G2M.Score") +
  ggtitle("Cell Cycle Scores")

print(p5)

# Correlation heatmap of QC metrics
qc_metrics <- sc_569@meta.data %>%
  select(nCount_SCT, nFeature_SCT, mitoRatio, , S.Score, G2M.Score)

qc_cor <- cor(qc_metrics, use = "complete.obs")
pheatmap(qc_cor, 
         display_numbers = TRUE, 
         number_format = "%.2f", 
         main = "Correlation between QC Metrics",
         color = colorRampPalette(c("blue", "white", "red"))(100))

dev.off()
