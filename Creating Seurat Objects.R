# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

process_gse_datasets <- function(gse_ids, base_dir) {
  # Initialize list for Seurat objects
  samples_1 <- list()
  
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
    
    if (length(matrix_files) == 0) {
      warning(sprintf("No matrix files found for %s", gse_id))
      next
    }
    
    cat(sprintf("Found %d matrix files\n", length(matrix_files)))
    
    # Process each matrix file
    for (matrix_file in matrix_files) {
      # Extract sample identifier from the matrix file
      file_prefix <- sub("_matrix\\.mtx(\\.gz)?$", "", basename(matrix_file))
      
      # Find the corresponding feature and barcode files
      feature_file <- all_files[grepl(paste0(file_prefix, "_features(\\.tsv|\\.txt)(\\.gz)?$"), all_files)]
      barcode_file <- all_files[grepl(paste0(file_prefix, "_barcodes(\\.tsv|\\.txt)(\\.gz)?$"), all_files)]
      
      # Ensure only one feature/barcode file per matrix
      if (length(feature_file) > 0) feature_file <- feature_file[1]
      if (length(barcode_file) > 0) barcode_file <- barcode_file[1]
      
      # Check if we have all required files
      if (length(feature_file) == 0 || length(barcode_file) == 0) {
        warning(sprintf("Missing feature or barcode file for matrix: %s", basename(matrix_file)))
        next
      }
      
      cat(sprintf("\nProcessing sample: %s\n", file_prefix))
      cat(sprintf("  Matrix: %s\n  Features: %s\n  Barcodes: %s\n", 
                  basename(matrix_file), basename(feature_file), basename(barcode_file)))
      
      # Create Seurat object
      tryCatch({
        counts <- ReadMtx(
          mtx = matrix_file,
          features = feature_file,
          cells = barcode_file,
          feature.column = 1
        )
        
        seurat_obj <- CreateSeuratObject(counts = counts, project = file_prefix)
        samples_1[[file_prefix]] <- seurat_obj
        cat(sprintf("  Created Seurat object for sample: %s\n", file_prefix))
      }, error = function(e) {
        warning(sprintf("Error creating Seurat object for %s: %s", file_prefix, e$message))
      })
    }
  }
  
  # Save all samples to an RDS file
  rds_path <- file.path(base_dir, "samples_1.rds")
  saveRDS(samples_1, file = rds_path)
  cat(sprintf("\nSaved all Seurat objects to: %s\n", rds_path))
  
  return(samples_1)
}

# Example usage
gse_ids <- c("GSE223128")  # Replace with your actual GSE IDs
base_dir <- "H:/ADB/mitoDynamics/single cell"  # Replace with your actual base directory

# Process datasets and get all Seurat objects
samples_1 <- process_gse_datasets(gse_ids, base_dir)

