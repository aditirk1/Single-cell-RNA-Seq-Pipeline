# RAW.tar Files with barcode, matrix and feature files for automated suerat list generation
# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(irGSEA)
library(scales)
library(patchwork)
library(SeuratData)
library(SeuratWrappers)
library(ComplexHeatmap)
library(AnnotationDbi)
library(future)
library(RColorBrewer)
#set working directory
setwd("/home/saini_lab/Documents/mitoDynamics") 
#annotations
ah <- AnnotationHub()
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
#Check versions of databases available
ahDb %>% 
  mcols()
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
#Download the appropriate Ensembldb database
edb <- ah[[id]]
# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

#View(annotations)                    
# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
# Explore biotypes
annotations$gene_biotype %>%
  factor() %>%
  levels()

# Extract IDs for mitochondrial genes
mt <- annotations %>% 
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_name)

process_gse_datasets <- function(base_dir) {
  # Find all *RAW.tar files in the base directory
  tar_files <- list.files(base_dir, pattern = ".*RAW\\.tar$", full.names = TRUE)
  # tar_files <- "GSE273937_RAW.tar"
  if (length(tar_files) == 0) {
    cat("No *RAW.tar files found in the specified directory.\n")
    return()
  }
  
  cat(sprintf("Found %d *RAW.tar files to process\n", length(tar_files)))
  
  for (tar_file in tar_files) {
    # Extract GSE ID from the tar file name
    gse_id <- sub("_RAW\\.tar$", "", basename(tar_file))
    cat(sprintf("\nProcessing file: %s (GSE ID: %s)\n", basename(tar_file), gse_id))

    # Create GSE directory if it doesn't exist
    gse_dir <- file.path(base_dir, gse_id)
    dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)

    # Create a temporary extraction directory
    extract_dir <- file.path(gse_dir, "temp_extract")
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

    # Extract tar file
    cat(sprintf("Extracting: %s\n", basename(tar_file)))
    untar(tar_file, exdir = extract_dir)

    # Copy the original tar file to the GSE directory if it's not already there
    if (!file.exists(file.path(gse_dir, basename(tar_file)))) {
      file.copy(tar_file, gse_dir)
      cat(sprintf("Copied tar file to: %s\n", file.path(gse_dir, basename(tar_file))))
    }

    # Find all relevant files in the extracted directory
    all_files <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)

    # Extract matrix files
    matrix_files <- all_files[grepl("matrix\\.mtx(\\.gz)?$", all_files)]
    
    cat(sprintf("Found %d matrix files\n", length(matrix_files)))
    
    # Initialize list to store all Seurat objects for this GSE ID
    seurat_list <- list()
    
    # Process each matrix file
    for (matrix_file in matrix_files) {
      # Extract sample identifier from the matrix file
      file_prefix <- sub("matrix\\.mtx(\\.gz)?$", "", basename(matrix_file))
      
      # Find the corresponding feature and barcode files
      feature_file <- all_files[grepl(paste0(file_prefix, "(features|genes)(\\.tsv|\\.txt)(\\.gz)?$"), all_files)][1]
      barcode_file <- all_files[grepl(paste0(file_prefix, "barcodes(\\.tsv|\\.txt)(\\.gz)?$"), all_files)][1]

      # meta_data_file <- all_files[grepl(paste0(file_prefix, "_bc_samples\\.csv(\\.gz)?$"), all_files)][1]
      # Check if required files exist
      if (is.na(feature_file) || is.na(barcode_file)) {
        cat(sprintf("Warning: Missing required files for sample %s, skipping\n", file_prefix))
        next
      }
      cat(sprintf("\nProcessing sample: %s\n", file_prefix))
      cat(sprintf("  Matrix: %s\n  Features: %s\n  Barcodes: %s\n", 
                  basename(matrix_file), basename(feature_file), basename(barcode_file)))
      
      # Read the count matrix
      counts <- Matrix::readMM(if(grepl("\\.gz$", matrix_file)) gzfile(matrix_file) else matrix_file)
      counts <- as(counts, "CsparseMatrix")
      
      # Read the features file
      genes <- readr::read_tsv(if(grepl("\\.gz$", feature_file)) gzfile(feature_file) else feature_file, 
                               col_names = FALSE, show_col_types = FALSE)
      gene_ids <- genes$X1
      
      # Read the barcodes file
      cell_ids <- readr::read_tsv(if(grepl("\\.gz$", barcode_file)) gzfile(barcode_file) else barcode_file, 
                                  col_names = FALSE, show_col_types = FALSE)$X1
      # Read the count matrix
      counts <- Matrix::readMM(if(grepl("\\.gz$", matrix_file)) gzfile(matrix_file) else matrix_file)
      if (is.null(counts) || length(counts) == 0) {
        stop(sprintf("Error: matrix file %s appears empty or failed to load", matrix_file))
      }
      counts <- as(counts, "CsparseMatrix")
      
      # Read features file
      genes <- readr::read_tsv(if(grepl("\\.gz$", feature_file)) gzfile(feature_file) else feature_file, 
                               col_names = FALSE, show_col_types = FALSE)
      gene_ids <- genes$X1
      
      # Read barcodes file
      cell_ids <- readr::read_tsv(if(grepl("\\.gz$", barcode_file)) gzfile(barcode_file) else barcode_file, 
                                  col_names = FALSE, show_col_types = FALSE)$X1
      
      # Check dimensions
      if (!inherits(counts, "Matrix")) {
        stop(sprintf("Error: matrix file %s did not load as a sparse matrix", matrix_file))
      }
      if (nrow(counts) != length(gene_ids)) {
        cat(sprintf("Warning: matrix row count (%d) != feature count (%d), trimming feature list\n", 
                    nrow(counts), length(gene_ids)))
        gene_ids <- gene_ids[1:nrow(counts)]
      }
      if (ncol(counts) != length(cell_ids)) {
        cat(sprintf("Warning: matrix column count (%d) != barcode count (%d), trimming barcode list\n", 
                    ncol(counts), length(cell_ids)))
        cell_ids <- cell_ids[1:ncol(counts)]
      }
      
      # Process metadata if available
      if (!is.na(meta_data_file)) {
        # Read metadata
        suppressMessages({
          map <- readr::read_csv(if(grepl("\\.gz$", meta_data_file)) gzfile(meta_data_file) else meta_data_file, 
                                 col_names = TRUE, show_col_types = FALSE)
        })
        colnames(map) <- c("sample_map", "sample")
        
        # Prepare metadata
        cell_ids3 <- sapply(strsplit(cell_ids, "_"), "[", 3)
        metadata <- data.frame(row.names = cell_ids, cells = cell_ids, sample_map = cell_ids3, stringsAsFactors = FALSE)
        metadata$sample <- map$sample[match(metadata$sample_map, map$sample_map)]
        
        # Try to map gene symbols if edb is available
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
      seurat_obj$mitopercentage <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
      
      # Standardize column names
      seurat_obj@meta.data <- dplyr::rename(seurat_obj@meta.data, 
                                            seq_folder = orig.ident,
                                            nUMI = nCount_RNA,
                                            nGene = nFeature_RNA)
      
      # Add Seurat object to the list with the sample name as the key
      seurat_list[[file_prefix]] <- seurat_obj
      cat(sprintf("  Added sample %s to the list\n", file_prefix))
    }
    
    # Save the entire list of Seurat objects for this GSE ID as a single RDS file
    if (length(seurat_list) > 0) {
      rds_filename <- file.path(gse_dir, paste0(gse_id, "_seurat_objects.rds"))
      saveRDS(seurat_list, file = rds_filename)
      cat(sprintf("\nSaved list of %d Seurat objects to: %s\n", length(seurat_list), rds_filename))
    } else {
      cat(sprintf("\nWarning: No Seurat objects were created for %s\n", gse_id))
    }
    
    # Clean up temporary extraction directory
    unlink(extract_dir, recursive = TRUE)
    cat(sprintf("Cleaned up temporary files for %s\n", gse_id))
    
    cat(sprintf("\nCompleted processing %s\n", gse_id))
  }
  
  cat("\nAll *RAW.tar files processed successfully.\n")
}

# Example usage
base_dir <- "/home/saini_lab/Documents/mitoDynamics/data"  # Replace with your actual base directory

# Process all *RAW.tar files and save consolidated Seurat object lists
process_gse_datasets(base_dir)
# Function to merge specific samples
merge_specific_samples <- function(samples) {
  # Merge sc5, sc6, and sc9
  sc_569 <- merge(x = samples[["GSM6940121_sc5"]], 
                  y = c(samples[["GSM6940122_sc6"]], samples[["GSM6940124_sc9"]]),
                  add.cell.id = c("tc1", "tc2", "tc3"))
  
  # Create metadata dataframe with group labels
  metadata <- sc_569@meta.data
  metadata$cells <- rownames(metadata)
  metadata$group <- NA
  metadata$group[which(stringr::str_detect(metadata$cells, "^tc1_"))] <- "tc1"
  metadata$group[which(stringr::str_detect(metadata$cells, "^tc2_"))] <- "tc2"
  metadata$group[which(stringr::str_detect(metadata$cells, "^tc3_"))] <- "tc3"
  
  sc_569@meta.data <- metadata
  sc_569@meta.data$sample[sc_569@meta.data$sample == "control"] <- "Control"
  valid_cells <- colnames(sc_569)[!is.na(sc_569@meta.data$sample)]
  sc_569 <- subset(sc_569, cells = valid_cells)
  # Clean up samples list
  samples[["sc_569"]] <- sc_569
  samples$GSM6940121_sc5 <- NULL
  samples$GSM6940122_sc6 <- NULL
  samples$GSM6940124_sc9 <- NULL

  return(samples)
}
# Merge the specific samples
samples <- merge_specific_samples(samples)
# Save all samples to an RDS file
path_sapmle <- file.path(base_dir, "samples.rds")

#sample specific naming 
names(samples)[names(samples) == "GSM6940120_sc2"] <- "sc_2"
names(samples)[names(samples) == "GSM6940123_sc7_8"] <- "sc_7_8"
