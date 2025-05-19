library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(gplots)
library(msigdbr)
library(pheatmap)
library(viridis)
library(ggrepel)

# Function to find pathway driver genes
FindPathwayDrivers <- function(seurat_obj, pathway_name, module_col, top_n = 10) {
  # Get correlation of each gene with the pathway score
  cat("Finding driver genes for", pathway_name, "pathway\n")
  
  # Extract expression data
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Get pathway scores
  pathway_scores <- seurat_obj@meta.data[[module_col]]
  
  # Calculate correlation for each gene
  gene_cors <- sapply(rownames(expr_data), function(gene) {
    tryCatch({
      cor(as.numeric(expr_data[gene,]), pathway_scores)
    }, error = function(e) NA)
  })
  
  # Get top correlated genes
  top_genes <- names(sort(gene_cors, decreasing = TRUE))[1:top_n]
  top_cors <- gene_cors[top_genes]
  
  # Return as data frame
  data.frame(
    gene = top_genes,
    correlation = top_cors,
    pathway = pathway_name,
    stringsAsFactors = FALSE
  )
}

# Create output directory if it doesn't exist
if (!dir.exists("output_plots")) {
  dir.create("output_plots")
  cat("Created output_plots directory\n")
}

# Get gene sets from MSigDB for key pathways
# Define the specific pathways we want to analyze
cat("Fetching MSigDB gene sets for pathway analysis...\n")


# Process each sample
for (i in names(samples)) {
  # Add module scores for pathway analysis
  cat("Calculating pathway scores for", i, "\n")
  for (pathway_name in names(pathway_gene_sets)) {
    # Get genes that are present in the dataset
    genes_in_data <- intersect(pathway_gene_sets[[pathway_name]], rownames(samples[[i]]))
    
    if (length(genes_in_data) < 5) {
      cat("WARNING: Too few genes found for pathway", pathway_name, "in sample", i, "(", length(genes_in_data), "genes). Skipping.\n")
      next
    }
    cat("  - Adding module score for", pathway_name, "pathway (", length(genes_in_data), "genes)\n")
    samples[[i]] <- AddModuleScore(
      object = samples[[i]],
      features = list(genes_in_data),
      name = pathway_name,
      seed = 123
    )
  }
}

# Generate visualizations and analysis for each sample
for (i in names(samples)) {
  # Save main UMAP plot with clusters
  p1 <- DimPlot(samples[[i]], reduction = "umap.harmony", label = TRUE)
  ggsave(paste0("output_plots/", i, "_UMAP_clusters.pdf"), p1 + ggtitle(paste("UMAP for", i)), width = 10, height = 8)
  
  # Save Elbow plot
  p2 <- ElbowPlot(samples[[i]], ndims = 50)
  ggsave(paste0("output_plots/", i, "_ElbowPlot.pdf"), p2 + ggtitle(paste("Elbow plot for", i)), width = 8, height = 6)
  
  # Create pathway score visualizations on UMAP
  cat("Generating pathway score visualizations for", i, "\n")
  pathway_plots <- list()
  
  # Create feature plots for each pathway
  for (pathway_name in names(pathway_gene_sets)) {
    module_col <- paste0(pathway_name, "1")
    
    if (module_col %in% colnames(samples[[i]]@meta.data)) {
      # Create and save feature plot
      p <- FeaturePlot(samples[[i]], features = module_col, 
                       min.cutoff = "q10", max.cutoff = "q90", label = TRUE) +
        scale_colour_gradientn(colors = viridis(100)) +
        ggtitle(paste0(i, ": ", pathway_name, " pathway score"))
      
      ggsave(paste0("output_plots/", i, "_", pathway_name, "_score.pdf"), p, width = 10, height = 8)
      pathway_plots[[pathway_name]] <- p
    }
  }
  
  # Create heatmap of pathway scores across clusters
  cat("Generating pathway score heatmap for", i, "\n")
  
  # Get pathway scores from metadata and average by cluster
  meta_data <- samples[[i]]@meta.data
  pathway_cols <- grep("^(hypoxia|emt|glycolysis|oxidative_phosphorylation|reactive_oxygen_species|apoptosis|inflammatory_response|p53_pathway|mitochondrial)1$", 
                       colnames(meta_data), value = TRUE)
  
  if (length(pathway_cols) > 0) {
    # Create a data frame with cluster and pathway scores
    pathway_data <- meta_data %>%
      select(seurat_clusters, all_of(pathway_cols)) %>%
      group_by(seurat_clusters) %>%
      summarize_all(mean)
    
    # Prepare data for heatmap
    pathway_matrix <- as.matrix(pathway_data[, -1])
    rownames(pathway_matrix) <- paste0("Cluster ", pathway_data$seurat_clusters)
    colnames(pathway_matrix) <- gsub("1$", "", pathway_cols)
    
    # Create heatmap
    pheatmap(
      pathway_matrix,
      main = paste0(i, ": Pathway Scores by Cluster"),
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      angle_col = 45,
      cellwidth = 20,
      cellheight = 20,
      border_color = NA,
      filename = paste0("output_plots/", i, "_pathway_heatmap.pdf"),
      width = 10,
      height = 8
    )
  } else {
    cat("No pathway scores available for heatmap in sample", i, "\n")
  }
  
  # Find top genes correlated with each pathway (pathway drivers)
  cat("Finding pathway driver genes for", i, "\n")
  all_drivers <- data.frame()
  
  for (pathway_name in names(pathway_gene_sets)) {
    module_col <- paste0(pathway_name, "1")
    
    if (module_col %in% colnames(samples[[i]]@meta.data)) {
      # Find pathway drivers
      drivers <- FindPathwayDrivers(samples[[i]], pathway_name, module_col, top_n = 15)
      all_drivers <- rbind(all_drivers, drivers)
      
      # Create and save visualization of top drivers
      p <- ggplot(drivers, aes(x = reorder(gene, correlation), y = correlation)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(x = "Gene", y = "Correlation with pathway score", 
             title = paste0(i, ": Top genes correlated with ", pathway_name, " pathway")) +
        theme_minimal()
      
      ggsave(paste0("output_plots/", i, "_", pathway_name, "_drivers.pdf"), p, width = 10, height = 8)
    }
  }
  
  # Save list of all drivers
  if (nrow(all_drivers) > 0) {
    write.csv(all_drivers, paste0("output_plots/", i, "_all_pathway_drivers.csv"), row.names = FALSE)
  }
}

# Perform differential expression analysis between clusters
cat("Performing differential expression analysis between clusters\n")
for (i in names(samples)) {
  cat("Finding cluster markers for", i, "\n")
  
  # Find all markers for every cluster
  all_markers <- FindAllMarkers(samples[[i]], 
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
  
  # Write results to CSV
  write.csv(all_markers, paste0("output_plots/", i, "_cluster_markers.csv"), row.names = FALSE)
  
  # Create heatmap of top cluster markers
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  
  # Plot top markers in a heatmap
  if (nrow(top_markers) > 0) {
    pdf(paste0("output_plots/", i, "_marker_heatmap.pdf"), width = 12, height = 12)
    DoHeatmap(samples[[i]], features = top_markers$gene) + 
      scale_fill_gradientn(colors = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
    
    # Create dot plot of top markers
    pdf(paste0("output_plots/", i, "_marker_dotplot.pdf"), width = 12, height = 10)
    DotPlot(samples[[i]], features = unique(top_markers$gene), cols = c("lightgrey", "blue")) + 
      RotatedAxis()
    dev.off()
  }
}

# Perform final comparison of all samples (if multiple)
if (length(samples) > 1) {
  cat("Performing comparative analysis across samples\n")
  
  # Compare pathway scores across samples
  all_sample_scores <- data.frame()
  
  for (i in names(samples)) {
    # Get pathway scores
    meta_data <- samples[[i]]@meta.data
    pathway_cols <- grep("^(hypoxia|emt|glycolysis|oxidative_phosphorylation|reactive_oxygen_species|apoptosis|inflammatory_response|p53_pathway|mitochondrial)1$", 
                         colnames(meta_data), value = TRUE)
    
    if (length(pathway_cols) > 0) {
      sample_scores <- meta_data %>%
        select(all_of(pathway_cols)) %>%
        summarize_all(mean) %>%
        mutate(sample = i)
      
      all_sample_scores <- rbind(all_sample_scores, sample_scores)
    }
  }
  
  # Create comparative heatmap if more than one sample
  if (nrow(all_sample_scores) > 1) {
    # Prepare data for heatmap
    score_matrix <- as.matrix(all_sample_scores[, !colnames(all_sample_scores) %in% "sample"])
    rownames(score_matrix) <- all_sample_scores$sample
    colnames(score_matrix) <- gsub("1$", "", colnames(score_matrix))
    
    # Create comparative heatmap
    pheatmap(
      score_matrix,
      main = "Pathway Scores Comparison Across Samples",
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      angle_col = 45,
      cellwidth = 20,
      cellheight = 20,
      border_color = NA,
      filename = "output_plots/all_samples_pathway_comparison.pdf",
      width = 10,
      height = 8
    )
  }
}

cat("Analysis complete! Results saved in 'output_plots' directory\n")

