# Create empty list to store plots
plots <- list()
if (!dir.exists("qc_metrics_SCT")) {
  dir.create("qc_metrics_SCT")
}
# Define metrics to plot
metrics_to_plot <- c("nUMI", "nGene", "log10GenesPerUMI", "mitoRatio")

for (i in names(samples)) {
  meta_data <- samples[[i]]@meta.data
  sample_name <- i  # Ensure sample_name is defined as the name of the sample
  
  # Create a custom directory for each sample
  dir.create(file.path("qc_metrics_SCT", sample_name), showWarnings = FALSE)
  
  cat("Processing QC plots for", sample_name, "\n")
  
  # Save the plots to the custom directory
  pdf(file.path("qc_metrics_SCT", sample_name, paste0("before_filt_", sample_name, ".pdf")), width = 10, height = 8)
  
  # Violin + boxplot for basic QC metrics
  print(meta_data %>%
          pivot_longer(cols = all_of(metrics_to_plot),
                       names_to = "metric", values_to = "value") %>%
          ggplot(aes(x = sample_name, y = value, fill = sample_name)) +
          geom_violin(alpha = 0.4) +
          geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
          facet_wrap(~metric, scales = "free_y", ncol = 2) +
          theme_bw() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) +
          labs(title = paste("Distribution of QC Metrics -", sample_name),
               x = "", y = "Value"))
  
  # Violin + boxplot grouped by 'group' if present
  if ("group" %in% colnames(meta_data)) {
    print(meta_data %>%
            pivot_longer(cols = all_of(metrics_to_plot),
                         names_to = "metric", values_to = "value") %>%
            ggplot(aes(x = group, y = value, fill = group)) +
            geom_violin(alpha = 0.4) +
            geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
            facet_wrap(~metric, scales = "free_y", ncol = 2) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste("Distribution of QC Metrics by Group -", sample_name),
                 x = "Group", y = "Value"))
  }
  
  # Feature count relationship plot (genes vs UMIs)
  if ("group" %in% colnames(meta_data)) {
    print(ggplot(meta_data, aes(x = nUMI, y = nGene)) +
            geom_point(aes(color = group), alpha = 0.6) +
            scale_color_brewer(palette = "Set1") +
            geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
            labs(title = paste("Gene vs UMI Counts -", sample_name),
                 x = "UMI Count", y = "Gene Count", color = "Group") +
            theme_bw() +
            scale_x_log10() + 
            scale_y_log10())
  } else {
    print(ggplot(meta_data, aes(x = nUMI, y = nGene)) +
            geom_point(aes(color = mitoRatio), alpha = 0.6) +
            scale_color_viridis_c() +
            geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
            labs(title = paste("Gene vs UMI Counts -", sample_name),
                 x = "UMI Count", y = "Gene Count", color = "Mito Ratio") +
            theme_bw() +
            scale_x_log10() + 
            scale_y_log10())
  }
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          scale_x_log10() + 
          theme_classic() +
          ylab("Cell density") +
          ggtitle(paste0("UMI Distribution - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = nGene, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          ggtitle(paste0("Gene Distribution - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          ggtitle(paste0("Mitochondrial % - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = log10GenesPerUMI, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          ggtitle(paste0("Log10 Genes Per UMI - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black") +
          stat_smooth(method = lm, aes(group = 1)) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          facet_wrap(~sample) +
          ggtitle(paste0("Gene vs UMI - ", sample_name)))
  # Plot data distributions
  print(samples[[i]]@meta.data %>%
          ggplot(aes(x=sample, fill=sample)) +
          geom_bar() +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("Cells Per Sample"))

  dev.off()  # Close the PDF device
  cat("Saved QC plots for", sample_name, "\n")
}

# Store the original cell counts before filtering
cells_before <- sapply(samples, ncol)

# Now perform filtering on each sample with custom thresholds
samples[["sc_2"]] <- subset(x = samples[["sc_2"]], 
                            subset = (nUMI >= 500 & nUMI <= 10000) & 
                              (nGene >= 250) & 
                              (log10GenesPerUMI > 0.88) & 
                              (mitoRatio < 0.1))

samples[["sc_7_8"]] <- subset(x = samples[["sc_7_8"]], 
                              subset = (nUMI >= 500 & nUMI <= 10000) & 
                                (nGene >= 300) & 
                                (log10GenesPerUMI > 0.88) & 
                                (mitoRatio < 0.1))

samples[["sc_569"]] <- subset(x = samples[["sc_569"]], 
                              subset = (nUMI >= 500) & 
                                (nGene >= 250 & nGene <= 10000) & 
                                (log10GenesPerUMI > 0.85) & 
                                (mitoRatio < 0.1))

# Calculate cells after filtering
cells_after <- sapply(samples, ncol)

# Process each sample and create plots
for (i in names(samples)) {
  # Get sample name
  sample_name <- i
  meta_data <- samples[[i]]@meta.data
  # Calculate statistics
  cells_removed <- cells_before[i] - cells_after[i]
  percent_removed <- round((cells_removed / cells_before[i]) * 100, 2)
  
  # Print filtering summary
  cat(paste0("Sample: ", sample_name, "\n"))
  cat(paste0(" Cells before filtering: ", cells_before[i], "\n"))
  cat(paste0(" Cells after filtering: ", cells_after[i], "\n"))
  cat(paste0(" Cells removed: ", cells_removed, " (", percent_removed, "%)\n\n"))
  
  # Set PDF filename to save it in the respective sample folder
  pdf_filename <- file.path("qc_metrics_SCT", sample_name, paste0("after_filt_", sample_name, ".pdf"))
  pdf(pdf_filename, width = 10, height = 8)
  
  # Violin + boxplot for basic QC metrics
  print(meta_data %>%
          pivot_longer(cols = all_of(metrics_to_plot),
                       names_to = "metric", values_to = "value") %>%
          ggplot(aes(x = sample_name, y = value, fill = sample_name)) +
          geom_violin(alpha = 0.4) +
          geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
          facet_wrap(~metric, scales = "free_y", ncol = 2) +
          theme_bw() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) +
          labs(title = paste("Distribution of QC Metrics -", sample_name),
               x = "", y = "Value"))
  
  # Violin + boxplot grouped by 'group' if present
  if ("group" %in% colnames(meta_data)) {
    print(meta_data %>%
            pivot_longer(cols = all_of(metrics_to_plot),
                         names_to = "metric", values_to = "value") %>%
            ggplot(aes(x = group, y = value, fill = group)) +
            geom_violin(alpha = 0.4) +
            geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
            facet_wrap(~metric, scales = "free_y", ncol = 2) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste("Distribution of QC Metrics by Group -", sample_name),
                 x = "Group", y = "Value"))
  }
  
  # Feature count relationship plot (genes vs UMIs)
  if ("group" %in% colnames(meta_data)) {
    print(ggplot(meta_data, aes(x = nUMI, y = nGene)) +
            geom_point(aes(color = group), alpha = 0.6) +
            scale_color_brewer(palette = "Set1") +
            geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
            labs(title = paste("Gene vs UMI Counts -", sample_name),
                 x = "UMI Count", y = "Gene Count", color = "Group") +
            theme_bw() +
            scale_x_log10() + 
            scale_y_log10())
  } else {
    print(ggplot(meta_data, aes(x = nUMI, y = nGene)) +
            geom_point(aes(color = mitoRatio), alpha = 0.6) +
            scale_color_viridis_c() +
            geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
            labs(title = paste("Gene vs UMI Counts -", sample_name),
                 x = "UMI Count", y = "Gene Count", color = "Mito Ratio") +
            theme_bw() +
            scale_x_log10() + 
            scale_y_log10())
  }
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          scale_x_log10() + 
          theme_classic() +
          ylab("Cell density") +
          ggtitle(paste0("UMI Distribution - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = nGene, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          ggtitle(paste0("Gene Distribution - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          ggtitle(paste0("Mitochondrial % - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(color = sample, x = log10GenesPerUMI, fill = sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          ggtitle(paste0("Log10 Genes Per UMI - ", sample_name)))
  
  print(samples[[i]]@meta.data %>% 
          ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black") +
          stat_smooth(method = lm, aes(group = 1)) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          facet_wrap(~sample) +
          ggtitle(paste0("Gene vs UMI - ", sample_name)))
  
  # Plot filtered data distributions
  print(samples[[i]]@meta.data %>%
    ggplot(aes(x=sample, fill=sample)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("Cells Per Sample"))

  dev.off() # Close the PDF device
  cat("Saved after filtering plots for", sample_name, "\n")
}
saveRDS(samples, file = "/home/saini_lab/Documents/mitoDynamics/Processed files/GSE223128/samples_pre-sct_filtered")
