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

process_gse_datasets <- function(gse_ids, base_dir) {
  # Initialize list for Seurat objects
  samples <- list()
  
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
      
      samples[[file_prefix]] <- seurat_obj
      cat(sprintf("  Created Seurat object for sample: %s\n", file_prefix))
    }
  }
  return(samples)
}

# Example usage
gse_ids <- c("GSE223128")  # Replace with your actual GSE IDs
base_dir <- "/home/saini_lab/Documents/mitoDynamics/data"  # Replace with your actual base directory

# Process datasets and get all Seurat objects
samples <- process_gse_datasets(gse_ids, base_dir)
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
