# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

process_gse_datasets <- function(gse_ids, base_dir) {
  # Initialize list for Seurat objects
  samples_2 <- list()
  
  for (gse_id in gse_ids) {
    cat(sprintf("\nProcessing GSE ID: %s\n", gse_id))
    
    # Create GSE directory if it doesn't exist
    gse_dir <- file.path(base_dir, paste0(gse_id, "_RAW"))
    dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Extract tar files if needed
    tar_files <- list.files(base_dir, pattern = paste0(gse_id, ".*\\.tar$"), full.names = TRUE)
    if (length(tar_files) > 0 && length(list.files(gse_dir)) == 0) {
      for (tar_file in tar_files) {
        untar(tar_file, exdir = gse_dir)
        cat(sprintf("Extracted: %s\n", basename(tar_file)))
      }
    }
    
    # Find all relevant files
    all_files <- list.files(gse_dir, full.names = TRUE, recursive = TRUE)
    # Extract matrix files
    matrix_files <- all_files[grepl("matrix\\.mtx(\\.gz)?$", all_files)]
    
    cat(sprintf("Found %d matrix files\n", length(matrix_files)))
    
    # Process each matrix file
    for (matrix_file in matrix_files) {
      # Extract sample identifier from the matrix file
      file_prefix <- sub("_matrix\\.mtx(\\.gz)?$", "", basename(matrix_file))
      # Find the corresponding feature and barcode files
      feature_file <- all_files[grepl(paste0(file_prefix, "_features(\\.tsv|\\.txt)(\\.gz)?$"), all_files)][1]
      barcode_file <- all_files[grepl(paste0(file_prefix, "_barcodes(\\.tsv|\\.txt)(\\.gz)?$"), all_files)][1]
      meta_data_file <- all_files[grepl(paste0(file_prefix, "_bc_samples\\.csv(\\.gz)?$"), all_files)][1]
      
      cat(sprintf("\nProcessing sample: %s\n", file_prefix))
      cat(sprintf("  Matrix: %s\n  Features: %s\n  Barcodes: %s\n", 
                  basename(matrix_file), basename(feature_file), basename(barcode_file)))
      
      # Read the count matrix
      counts <- Matrix::readMM(if(grepl("\\.gz$", matrix_file)) gzfile(matrix_file) else matrix_file)
      counts <- as(counts, "CsparseMatrix")
      
      # Read the features file
      genes <- readr::read_tsv(if(grepl("\\.gz$", feature_file)) gzfile(feature_file) else feature_file, col_names = FALSE)
      gene_ids <- genes$X1
      
      # Read the barcodes file
      cell_ids <- readr::read_tsv(if(grepl("\\.gz$", barcode_file)) gzfile(barcode_file) else barcode_file, col_names = FALSE)$X1
      
      # Process metadata if available
      if (!is.na(meta_data_file)) {
        # Read metadata
        suppressMessages({
          map <- readr::read_csv(if(grepl("\\.gz$", meta_data_file)) gzfile(meta_data_file) else meta_data_file, col_names = TRUE, show_col_types = FALSE)
        })
        colnames(map) <- c("sample_map", "sample")
        
        # Prepare metadata
        cell_ids3 <- sapply(strsplit(cell_ids, "_"), "[", 3)
        metadata <- data.frame(row.names = cell_ids, cells = cell_ids, sample_map = cell_ids3, stringsAsFactors = FALSE)
        metadata$sample <- map$sample[match(metadata$sample_map, map$sample_map)]
        
        # Try to map gene symbols if edb is available - FIX FOR ENSDB SELECT ISSUE
        if (exists("edb", envir = .GlobalEnv)) {
          # Use ensembldb's select method instead of dplyr's
          ann <- ensembldb::select(get("edb", envir = .GlobalEnv), keys = gene_ids, keytype = 'GENEID', columns = "SYMBOL")
          gene_id <- ifelse((gene_ids %in% ann$GENEID), ifelse(ann$SYMBOL == "", ann$GENEID, ann$SYMBOL), gene_ids)
        } else {
          gene_id <- gene_ids
        }
        gene_id <- make.unique(gene_id, sep = ".")
        
        # Set row and column names
        rownames(counts) <- gene_id
        colnames(counts) <- cell_ids
        
        # Create Seurat object with metadata
        seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata,
                                         min.features = 100, project = file_prefix)
      } else {
        # If no metadata, use simpler approach
        rownames(counts) <- make.unique(gene_ids, sep = ".")
        colnames(counts) <- cell_ids
        
        seurat_obj <- CreateSeuratObject(counts = counts, project = file_prefix)
      }
      
      # Add additional metrics
      seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
      seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
      seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
      
      samples_2[[file_prefix]] <- seurat_obj
      cat(sprintf("  Created Seurat object for sample: %s\n", file_prefix))
    }
    
    # Perform merging for specific samples if they exist (sc5, sc6, sc9)
    sample_pattern <- gsub("^GSE", "GSM", gse_id)
    sc5_key <- grep(paste0(sample_pattern, ".*_sc5$"), names(samples_2), value = TRUE)
    sc6_key <- grep(paste0(sample_pattern, ".*_sc6$"), names(samples_2), value = TRUE)
    sc9_key <- grep(paste0(sample_pattern, ".*_sc9$"), names(samples_2), value = TRUE)
    
    if (length(sc5_key) > 0 && length(sc6_key) > 0 && length(sc9_key) > 0) {
      cat("\nMerging sc5, sc6, and sc9 samples...\n")
      
      # Merge sc5, sc6, and sc9
      merged_seurat <- merge(x = samples_2[[sc5_key]], 
                             y = c(samples_2[[sc6_key]], samples_2[[sc9_key]]),
                             add.cell.id = c("tc1", "tc2", "tc3"))
      
      # Update metrics for merged object
      merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
      merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
      merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
      
      # Create metadata dataframe with group labels
      metadata <- merged_seurat@meta.data
      metadata$cells <- rownames(metadata)
      metadata$group <- NA
      metadata$group[which(stringr::str_detect(metadata$cells, "^tc1_"))] <- "tc1"
      metadata$group[which(stringr::str_detect(metadata$cells, "^tc2_"))] <- "tc2"
      metadata$group[which(stringr::str_detect(metadata$cells, "^tc3_"))] <- "tc3"
      
      # Rename columns for consistency
      metadata <- dplyr::rename(metadata, 
                                seq_folder = orig.ident,
                                nUMI = nCount_RNA,
                                nGene = nFeature_RNA)
      
      # Update metadata in merged object
      merged_seurat@meta.data <- metadata
      
      # Add to samples_2 list
      samples_2[["merged_sc5_sc6_sc9"]] <- merged_seurat
      cat("  Created merged Seurat object (sc5, sc6, sc9)\n")
    }
  }
  
  # Save all samples to an RDS file
  rds_path <- file.path(base_dir, "samples_2.rds")
  saveRDS(samples_2, file = rds_path)
  cat(sprintf("\nSaved all Seurat objects to: %s\n", rds_path))
  
  return(samples_2)
}
# Example usage
gse_ids <- c("GSE223128")  # Replace with your actual GSE IDs
base_dir <- "H:/ADB/mitoDynamics/single cell"  # Replace with your actual base directory

# Process datasets and get all Seurat objects
samples_2 <- process_gse_datasets(gse_ids, base_dir)
