## SCT NORMALIZATION
# 1️⃣ Set up parallel processing with memory limits
plan(strategy = "multicore", workers = 2)  # Adjust workers to your CPU cores
options(future.globals.maxSize = 100 * 1024^3)  # ~100 GB (adjust if needed)

# 2️⃣ Process each sample with SCTransform + PCA
for (i in names(samples)) {
  cat("Processing", i, "\n")
  
  # Clear memory aggressively before each sample
  gc(verbose = FALSE, full = TRUE)
  
  samples[[i]] <- SCTransform(samples[[i]], assay = "RNA", 
                              vars.to.regress = "G2M.Score", 
                              verbose = FALSE,)
  
  hvgs <- VariableFeatures(samples[[i]])
  
  samples[[i]] <- RunPCA(samples[[i]], assay = "SCT", 
                         features = hvgs, npcs = 30, verbose = FALSE)
}

# =======================================
# 3️⃣ Complete the processing (UMAP + clustering)
# =======================================
for (i in names(samples)) {
  cat("Completing processing for", i, "with UMAP and clustering\n")
  
  npcs <- 30
  
  samples[[i]] <- RunUMAP(samples[[i]], dims = 1:npcs, assay = "SCT", verbose = FALSE)
  samples[[i]] <- FindNeighbors(samples[[i]], dims = 1:npcs, assay = "SCT", verbose = FALSE)
  samples[["sc_569"]] <- FindClusters(samples[["sc_569"]], resolution = 0.4, algorithm = 4, verbose = FALSE)
  
  # Save each fully processed object to disk to avoid memory pileup
  saveRDS(samples[[i]], file = paste0("fully_processed_", i, ".rds"))

  gc(verbose = FALSE, full = TRUE)
}

#=======================================
# 2. Memory-efficient visualization functions
#=======================================
visualize_sample <- function(sample, sample_name, output_dir = "plots") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Base filename
  base_filename <- file.path(output_dir, sample_name)
  
  # 1. UMAP plot with clusters (no rasterization)
  p1 <- DimPlot(sample, reduction = "umap", group.by = "group", 
                label = TRUE, raster = FALSE) +
    ggtitle(paste0(sample_name, " - Clusters"))
  ggsave(paste0(base_filename, "_clusters.pdf"), p1, width = 8, height = 7, dpi = 300, device = "pdf")
  
  # 2. QC metrics
  if ("nCount_SCT" %in% colnames(sample@meta.data) && 
      "nFeature_SCT" %in% colnames(sample@meta.data)) {
    
    p2 <- VlnPlot(sample, features = c("nFeature_SCT", "nCount_SCT"), 
                  ncol = 2, pt.size = 0) +
      ggtitle(paste0(sample_name, " - QC Metrics"))
    ggsave(paste0(base_filename, "_qc.pdf"), p2, width = 10, height = 5, dpi = 300, device = "pdf")
    
    p3 <- FeatureScatter(sample, feature1 = "nCount_SCT", feature2 = "nFeature_SCT") +
      ggtitle(paste0(sample_name, " - Feature-Count Relationship"))
    ggsave(paste0(base_filename, "_feature_count.pdf"), p3, width = 7, height = 6, dpi = 300, device = "pdf")
  }
  
  # 3. Elbow plot
  p4 <- ElbowPlot(sample, ndims = 30) +
    ggtitle(paste0(sample_name, " - PCA Elbow Plot"))
  ggsave(paste0(base_filename, "_elbow.pdf"), p4, width = 7, height = 6, dpi = 300, device = "pdf")
  
  return(list(
    clusters = paste0(base_filename, "_clusters.pdf"),
    qc = paste0(base_filename, "_qc.pdf"),
    feature_count = paste0(base_filename, "_feature_count.pdf"),
    elbow = paste0(base_filename, "_elbow.pdf")
  ))
}

for(i in names(samples)){
  visualize_sample(samples[[i]], i)
}

# Function to visualize marker genes for a sample
visualize_markers <- function(sample, sample_name, output_dir = "plots", 
                              n_markers_per_cluster = 50) {
  
  # Create marker directory
  marker_dir <- file.path(output_dir, "markers")
  if (!dir.exists(marker_dir)) {
    dir.create(marker_dir, recursive = TRUE)
  }
  
  # Find markers for each cluster
  cat("Finding markers for", sample_name, "\n")
  
  # For memory efficiency, find markers for one cluster at a time
  clusters <- unique(sample$harmony_clusters_res0.2)
  all_markers <- data.frame()
  
  for (cluster in clusters) {
    gc(verbose = FALSE, full = TRUE)  # Clear memory before each cluster
    
    # Find markers for this cluster only
    cluster_markers <- FindMarkers(sample, 
                                   ident.1 = cluster, 
                                   min.pct = 0.25,
                                   only.pos = TRUE,
                                   max.cells.per.ident = 500)  # Limit cells for efficiency
    
    if (nrow(cluster_markers) > 0) {
      # Add cluster info
      cluster_markers$cluster <- cluster
      cluster_markers$gene <- rownames(cluster_markers)
      
      # Store the top markers
      top_markers <- cluster_markers %>%
        arrange(p_val_adj) %>%
        head(n_markers_per_cluster)
      
      all_markers <- rbind(all_markers, top_markers)
    }
  }
  
  # Save markers table
  marker_file <- file.path(marker_dir, paste0(sample_name, "_markers.csv"))
  write.csv(all_markers, marker_file)
  
  # Plot top markers per cluster
  for (cluster in unique(all_markers$cluster)) {
    gc(verbose = FALSE, full = TRUE)  # Clear memory
    
    cluster_genes <- all_markers %>% 
      filter(cluster == !!cluster) %>% 
      pull(gene)
    
    if (length(cluster_genes) > 0) {
      # Create feature plot for these genes
      p <- FeaturePlot(sample, 
                       features = cluster_genes[1:min(length(cluster_genes), 4)],
                       ncol = 2)
      
      ggsave(file.path(marker_dir, paste0(sample_name, "_cluster", cluster, "_markers.png")), 
             p, width = 10, height = 8, dpi = 300)
      
      # Create violin plot
      v <- VlnPlot(sample, 
                   features = cluster_genes[1:min(length(cluster_genes), 2)],
                   ncol = 2, pt.size = 0)
      
      ggsave(file.path(marker_dir, paste0(sample_name, "_cluster", cluster, "_vlnplot.png")), 
             v, width = 10, height = 5, dpi = 300)
    }
  }
  
  return(marker_file)
}

#=======================================
# 3. Generate visualizations efficiently
#=======================================
# Process each sample for visualization, one at a time
for (i in names(samples)) {
  cat("Visualizing", i, "\n")
  
  # Clear memory before visualization
  gc(verbose = FALSE, full = TRUE)
  
  # # Basic visualizations
  # visualize_sample(samples[[i]], i)
  
  # Marker gene visualizations (more memory intensive)
  visualize_markers(samples[[i]], i, n_markers_per_cluster = 2)
  
  gc(verbose = FALSE, full = TRUE)
  
  cat("Completed visualization for", i, "\n\n")
}

#### ABOVE HAS BEEN RUN SUCCESSFULLY.
# -------Pathways---------
# Get Hallmark gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(f = .$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# Get key pathway gene sets from HALLMARK collection
pathway_gene_sets <- list(
  hypoxia = hallmark_gene_sets$HALLMARK_HYPOXIA,
  emt = hallmark_gene_sets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  glycolysis = hallmark_gene_sets$HALLMARK_GLYCOLYSIS,
  oxidative_phosphorylation = hallmark_gene_sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
  reactive_oxygen_species = hallmark_gene_sets$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  apoptosis = hallmark_gene_sets$HALLMARK_APOPTOSIS,
  inflammatory_response = hallmark_gene_sets$HALLMARK_INFLAMMATORY_RESPONSE,
  p53_pathway = hallmark_gene_sets$HALLMARK_P53_PATHWAY
)

# Use existing MitoCarta pathways
#load mitocarta
cat("Adding MitoCarta mitochondrial pathways\n")

# Process mitochondrial pathways from MitoCarta data
mito_pathway_names <- unique(mito_pathways$MitoPathway)

# Add each mitochondrial pathway separately
for (name in mito_pathway_names) {
  # Skip if empty pathway name
  if (is.na(name) || name == "" || is.null(name)) {
    next
  }
  # Get genes for this pathway
  gene_strings <- mito_pathways %>%
    dplyr::filter(MitoPathway == name) %>%
    pull(Genes)
  # Split comma-separated gene strings and clean
  mito_genes <- unlist(lapply(gene_strings, function(x) {
    # Handle potential NAs
    if (is.na(x) || x == "") return(character(0))
    # Split by comma and trim whitespace
    genes <- unlist(strsplit(x, ","))
    trimws(genes)
  }))
  
  # Remove duplicates
  mito_genes <- unique(mito_genes)
  
  # Skip if too few genes
  if (length(mito_genes) < 5) {
    cat("Skipping mitochondrial pathway", name, "- too few genes\n")
    next
  }
  
  # Create clean pathway name for use in code
  clean_path_name <- gsub("[^[:alnum:]]", "_", tolower(name))
  clean_path_name <- paste0("mito_", clean_path_name)
  
  # Add to pathway gene sets
  pathway_gene_sets[[clean_path_name]] <- mito_genes
  cat("  - Added", name, "with", length(mito_genes), "genes\n")
}
