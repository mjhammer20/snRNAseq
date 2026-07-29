# Load Libraries
library(MAGMA.Celltyping)
library(Matrix)
library(scKirby)
library(anndataR)
library(foreach)
library(doParallel)
library(biomaRt)
library(argparse)

# Function to temporarily change working directory for an expression and then revert back
with_workdir <- function(dir, expr) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(dir)
  eval.parent(substitute(expr))
}

# Function to configure logging to file
  log_line <- function(msg) {
    cat(msg, "\n")
    flush.console()
  }


# Function to transform GWAS summary statistics to the format required by MAGMA, including annotation with rsIDs using biomaRt in parallel chunks
transform_gwas_sumstats <- function(gwas_file, output_file, chunk_size = 50000, n_workers = 4) {
  
  # Load the GWAS summary statistics
  gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t")

  # Annotate gwas sum stats with rsIDs using biomaRt (chunked + parallel)
  log_line("Annotating GWAS summary statistics with rsIDs using biomaRt (chunked)\n")
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

  total_rows <- nrow(gwas_data)
  chunk_ids <- split(seq_len(total_rows), ceiling(seq_len(total_rows) / chunk_size))
  log_line(sprintf("Total rows: %d | Chunks: %d | Chunk size: %d\n", total_rows, length(chunk_ids), chunk_size))

  # Setup parallel backend
  cl <- makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  snp_annotations_list <- foreach::foreach(i = seq_along(chunk_ids), .combine = rbind, .packages = "biomaRt") %dopar% {
    idx <- chunk_ids[[i]]
    chunk <- gwas_data[idx, ]

    # Progress logging (per worker)
    log_line(sprintf("  [Chunk %d/%d] Querying %d SNPs\n", i, length(chunk_ids), nrow(chunk)))

    getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),
      filters = c("chr_name", "start"),
      values = list(chunk$chromosome, chunk$base_pair_position),
      mart = ensembl
    )
  }

  # Merge the original GWAS data with the annotations to get rsIDs
  magma_gwas <- gwas_data %>%
    left_join(snp_annotations_list, by = c("chromosome" = "chr_name", "base_pair_position" = "start")) %>%
    rename(SNP = refsnp_id, CHR = chromosome, BP = base_pair_position, P = p_value, BETA = beta, SE = standard_error) %>%
    select(SNP, CHR, BP, P, BETA, SE)

  write.table(magma_gwas, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Function to prepare GWAS summary statistics for MAGMA, including annotation with rsIDs and mapping SNPs to genes
prepare_gwas_data <- function(gwas_summary_stats, magma_gwas_summary_stats, upstream_kb, downstream_kb) {
  log_line("Preparing GWAS summary statistics for MAGMA\n")
  
  # Transform the GWAS summary statistics to the format required by MAGMA
  log_line(sprintf("Input GWAS summary stats: %s\n", gwas_summary_stats))
  log_line(sprintf("Output MAGMA-formatted summary stats: %s\n", magma_gwas_summary_stats))
  transform_gwas_sumstats(gwas_summary_stats, magma_gwas_summary_stats)

  # Calculate total N
  cases <- 63555
  proxy_cases <- 17700
  controls <- 1746386
  total_n <- cases + proxy_cases + controls

  # Map SNPs to genes for MAGMA
  log_line("Mapping SNPs to genes for MAGMA\n")
  genes_out_path <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = magma_gwas_summary_stats,
    genome_build = "hg38",
    upstream_kb = upstream_kb,
    downstream_kb = downstream_kb,
    N = total_n
  )
  log_line(sprintf("Output gene mapping file: %s\n", genes_out_path))

  return(genes_out_path)
}

# Function to generate CTD from anndata
generate_ctd <- function(adata_prefix, region, cell_type_assignment_key) {
  # Load anndata
  log_line("Generating Cell Type Data (CTD) from anndata\n")
  adata_file <- file.path(paste0(adata_prefix, "_",region, ".h5ad"))
  log_line(sprintf("Reading anndata file: %s\n", adata_file))
  adata <- anndataR::read_h5ad(adata_file)

  # Debug: Check anndata structure
  log_line(sprintf("AnnData object info for region %s:\n", region))
  log_line(sprintf("  X shape (cells x genes): %s\n", paste(dim(adata$X), collapse = " x ")))
  log_line(sprintf("  Gene names (first 20): %s\n", paste(head(colnames(adata$X), 20), collapse = ", ")))

  # Debug: Check cell type annotations
  log_line(sprintf("  Unique cell types: %d\n", length(unique(adata$obs[[cell_type_assignment_key]]))))
  log_line(sprintf("  NA values in %s: %d\n", cell_type_assignment_key, sum(is.na(adata$obs[[cell_type_assignment_key]]))))
  log_line(sprintf("  Cell type distribution:\n"))
  cell_type_table <- table(adata$obs[[cell_type_assignment_key]])
  print(cell_type_table)

  # FILTER: Mark rare cell types as NA instead of subsetting
  min_cells_per_type <- 10
  rare_types <- names(cell_type_table[cell_type_table < min_cells_per_type])
  if (length(rare_types) > 0) {
    log_line(sprintf("  Filtering out %d rare cell types with < %d cells: %s\n", 
                length(rare_types), min_cells_per_type, paste(rare_types, collapse = ", ")))
    
    # Mark rare cell types as NA - keeps anndata structure intact
    adata$obs[[cell_type_assignment_key]][adata$obs[[cell_type_assignment_key]] %in% rare_types] <- NA
  }

  # Ensure X is cells x genes and sparse BEFORE any alignment checks
  log_line(sprintf("  X class: %s\n", paste(class(adata$X), collapse = ", ")))
  x_nrow <- tryCatch(nrow(adata$X), error = function(e) NA_integer_)
  x_ncol <- tryCatch(ncol(adata$X), error = function(e) NA_integer_)
  obs_nrow <- tryCatch(nrow(adata$obs), error = function(e) NA_integer_)
  log_line(sprintf("  X dim via nrow/ncol: %s x %s\n", x_nrow, x_ncol))
  log_line(sprintf("  obs rows: %s\n", obs_nrow))

  if (is.na(x_nrow) || is.na(x_ncol) || is.na(obs_nrow)) {
    stop("Unexpected shape for X (nrow/ncol not available)")
  }

  if (obs_nrow != x_nrow && obs_nrow == x_ncol) {
    log_line("  Transposing X to match obs rows\n")
    adata$X <- Matrix::t(adata$X)
    x_nrow <- nrow(adata$X)
  }
  if (obs_nrow != x_nrow) {
    stop(sprintf("Unexpected shape for X after alignment (obs rows=%s, X rows=%s)", obs_nrow, x_nrow))
  }

  tryCatch({
    adata$X <- as(adata$X, "CsparseMatrix")
  }, error = function(e) { stop(sprintf("FAIL: as(CsparseMatrix): %s", e$message)) })

  tryCatch({
    adata$X <- as(adata$X, "dgCMatrix")
  }, error = function(e) { stop(sprintf("FAIL: as(dgCMatrix): %s", e$message)) })

  # Ensure alignment and drop NA cells explicitly
  if (nrow(adata$obs) != nrow(adata$X)) {
    stop("Mismatch between obs rows and X rows")
  }

  # Identify cells with valid cell type annotations (non-NA) for subsetting
  tryCatch({
    keep_cells <- !is.na(adata$obs[[cell_type_assignment_key]])
  }, error = function(e) { stop(sprintf("FAIL: keep_cells: %s", e$message)) })
  
  if (length(keep_cells) != nrow(adata$X)) {
    stop(sprintf("keep_cells length (%s) != X rows (%s)", length(keep_cells), nrow(adata$X)))
  }

  # Grab cell indices with valid cell type annotations (non-NA) for subsetting
  keep_idx <- which(keep_cells)

  if (length(keep_idx) == 0) {
    stop(sprintf("No cells with valid %s remain for region %s", cell_type_assignment_key, region))
  }

  # Subset AnnData in one step to keep X/obs aligned
  tryCatch({
    adata <- adata[keep_idx, , drop = FALSE]
  }, error = function(e) { stop(sprintf("FAIL: adata subset: %s", e$message)) })

  # Convert view to concrete AnnData before modifying obs
  tryCatch({
    adata <- adata$as_InMemoryAnnData()
  }, error = function(e) { stop(sprintf("FAIL: as_InMemoryAnnData: %s", e$message)) })

  tryCatch({
    adata$obs[[cell_type_assignment_key]] <- droplevels(factor(adata$obs[[cell_type_assignment_key]]))
  }, error = function(e) { stop(sprintf("FAIL: droplevels: %s", e$message)) })

  log_line(sprintf("  Remaining cells with valid cell type: %d\n", nrow(adata$X)))
  log_line(sprintf("  Remaining cell types: %d\n", length(unique(adata$obs[[cell_type_assignment_key]]))))

  if (!cell_type_assignment_key %in% colnames(adata$obs)) {
    stop(sprintf("%s column missing", cell_type_assignment_key))
  }

  # Normalize obs/var names and annotations before CTD
  if (is.null(rownames(adata$obs))) {
    rownames(adata$obs) <- paste0("cell_", seq_len(nrow(adata$obs)))
  }
  if (is.null(colnames(adata$X))) {
    if (!is.null(adata$var_names)) {
      colnames(adata$X) <- as.character(adata$var_names)
    } else if (!is.null(adata$var$gene_ids)) {
      colnames(adata$X) <- as.character(adata$var$gene_ids)
    } else {
      stop("No gene names found for X columns")
    }
  }
  adata$obs[[cell_type_assignment_key]] <- as.character(adata$obs[[cell_type_assignment_key]])

  # Compatibility shim for orthogene::aggregate_rows signature
  ortho_ns <- asNamespace("orthogene")
  if (exists("aggregate_rows", envir = ortho_ns, inherits = FALSE)) {
    orig_fun <- get("aggregate_rows", envir = ortho_ns)
    if (!"as_delayedarray" %in% names(formals(orig_fun))) {
      unlockBinding("aggregate_rows", ortho_ns)
      assign(
        "aggregate_rows",
        function(..., as_delayedarray = NULL) orig_fun(...),
        envir = ortho_ns
      )
      lockBinding("aggregate_rows", ortho_ns)
    }
  } else {
    stop("orthogene::aggregate_rows not found in namespace")
  }

  # Create CTD from anndata
  log_line(sprintf("Creating CTD from anndata for region: %s\n", region))
  tryCatch({
    ctd <- scKirby::anndata_to_ctd(
      obj = adata, 
      annotLevels = list(adata$obs[[cell_type_assignment_key]]),
      input_species = "human",
      output_species = "human",
      dropNA = TRUE,
      force_standardise = TRUE
    )
    
    # Verify CTD has data
    if (is.null(ctd[[1]]$specificity) || nrow(ctd[[1]]$specificity) == 0) {
      stop(sprintf("CTD for region %s is empty - no specificity data\n", region))
    }

    log_line(sprintf("CTD created successfully with %d genes and %d cell types\n",
                nrow(ctd[[1]]$specificity), ncol(ctd[[1]]$specificity)))

    # Clean up
    rm(adata)
    gc()
    
    return(ctd)
  }, error = function(e) {
    log_line(sprintf("ERROR creating CTD: %s\n", e$message))
    stop(e)
  })
}

# Function to perform cell type enrichment analysis
perform_cell_type_enrichment <- function(gwas_base, ctd, region, output_dir, upstream_kb, downstream_kb, results_prefix) {
  log_line(sprintf("Performing cell type enrichment analysis for region: %s\n", region))

  # Define gwas sum stats path
  gwas_sumstats_path <- file.path(output_dir, gwas_base)

  # If ctd is a wrapper list, extract the actual CTD
  if (is.list(ctd) && !is.null(ctd$ctd)) {
    ctd <- ctd$ctd
  }

  # Run linear association test using MAGMA.Celltyping
  log_line(sprintf("Running linear association test for region: %s\n", region))
  enrichment_results_linear <- tryCatch({
    with_workdir(output_dir, {
      MAGMA.Celltyping::calculate_celltype_associations(
        ctd = ctd,
        ctd_levels = 1,
        gwas_sumstats_path = gwas_sumstats_path,
        ctd_species = "human",
        analysis_name = region,
        force_new = TRUE, 
        EnrichmentMode = "Linear",
        upstream_kb = upstream_kb,
        downstream_kb = downstream_kb
      )
    })
  }, error = function(e) {
    log_line(sprintf("ERROR in linear association test for region %s: %s\n", region, e$message))
    stop(e)
  })

  log_line("Linear association test complete\n")

  # Run top 10% association test using MAGMA.Celltyping
  log_line(sprintf("Running top 10%% association test for region: %s\n", region))
  enrichment_results_top10 <- tryCatch({
    with_workdir(output_dir, {
      MAGMA.Celltyping::calculate_celltype_associations(
        ctd = ctd,
        ctd_levels = 1,
        gwas_sumstats_path = gwas_sumstats_path,
        ctd_species = "human",
        analysis_name = region,
        force_new = TRUE, 
        EnrichmentMode = "Top 10%",
        upstream_kb = upstream_kb,
        downstream_kb = downstream_kb
      )
    })
  }, error = function(e) {
    log_line(sprintf("ERROR in top 10%% association test for region %s: %s\n", region, e$message))
    stop(e)
  })
  log_line("Top 10% association test complete\n")

  enrichment_results <- tryCatch({
    MAGMA.Celltyping::merge_magma_results(
      ctAssoc1 = enrichment_results_linear,
      ctAssoc2 = enrichment_results_top10
    )
  }, error = function(e) {
    log_line(sprintf("ERROR merging enrichment results for region %s: %s\n", region, e$message))
    stop(e)
  })

  # Save enrichment results to file
  results_fp = file.path(output_dir, paste0(results_prefix, "_", region, ".tsv"))
  log_line(sprintf("Saving enrichment results to %s", results_fp))
  write.table(enrichment_results, results_fp, sep = "\t", row.names = FALSE, quote = FALSE)

  return(enrichment_results)
}

# Function to plot results
plot_results <- function(enrichment_results, ctd, region, output_dir, gwas_base) {
  
  # Plot merged results
  log_line(sprintf("Plotting enrichment associations for region: %s\n", region))
  tryCatch({
    MAGMA.Celltyping::plot_celltype_associations(
      ctAssocs = enrichment_results,
      ctd = ctd,
      fileTag= paste0(region), 
      figsDir = file.path(output_dir, "MAGMA_Figures", gwas_base)
    )
  }, error = function(e) {
    log_line(sprintf("ERROR plotting enrichment associations for region %s: %s\n", region, e$message))
  })
  log_line("Plotting completed successfully\n")
}

main_precomputed <- function(ref_tax_output_dir, gwas_base, adata_prefix, regions, upstream_kb, downstream_kb, cell_type_assignment_key, results_prefix) {

  # Loop through regions and perform analysis
  results_list <- list()
  
  for (i in seq_along(regions)) {
    region <- regions[i]
    
    tryCatch({
      log_line(sprintf("[%d/%d] Processing region: %s\n", i, length(regions), region))
      log_line(sprintf("[%s] Processing region: %s\n", Sys.time(), region))

      # Generate CTD for the region
      ctd <- generate_ctd(adata_prefix, region, cell_type_assignment_key)
      log_line(sprintf("[%s] CTD generated\n", Sys.time()))
      
      # Perform cell type enrichment analysis
      enrichment_results <- perform_cell_type_enrichment(gwas_base, ctd, region, ref_tax_output_dir, upstream_kb, downstream_kb, results_prefix)
      log_line(sprintf("[%s] Enrichment analysis complete\n", Sys.time()))

      # Plot results
      plot_results(enrichment_results, ctd, region, ref_tax_output_dir, gwas_base)
      log_line(sprintf("[%s] Plotting complete\n", Sys.time()))

      # Save results
      saveRDS(enrichment_results, file.path(ref_tax_output_dir, paste0("precomputed_enrichment_results_", region, ".rds")))
      write.table(enrichment_results, file.path(ref_tax_output_dir, paste0("precomputed_enrichment_results_", region, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
      log_line(sprintf("[%s] Results saved\n", Sys.time()))
      
      # Clean up
      rm(ctd, enrichment_results)
      gc()
      
      results_list[[region]] <- "SUCCESS"
      
    }, error = function(e) {
      log_line(sprintf("[%s] ERROR: %s\n", Sys.time(), e$message))
      log_line(sprintf("ERROR processing region %s: %s\n", region, e$message))
      results_list[[region]] <<- paste("FAILED:", e$message)
    })
  }

  log_line("Analysis complete\n")
  log_line(sprintf("Successful regions: %d\n", sum(results_list == "SUCCESS")))
  log_line(sprintf("Failed regions: %d\n", sum(results_list != "SUCCESS")))
}

main <- function(ref_tax_output_dir, gwas_summary_stats, magma_gwas_summary_stats, adata_prefix, regions, upstream_kb, downstream_kb, cell_type_assignment_key, results_prefix) {

  # Prepare GWAS summary statistics for MAGMA
  genes_out_path <- prepare_gwas_data(gwas_summary_stats, magma_gwas_summary_stats, upstream_kb, downstream_kb)
  genes_out_path <- normalizePath(genes_out_path, mustWork = TRUE)

  # Use the formatted GWAS filename (no .genes.out)
  gwas_base <- basename(magma_gwas_summary_stats)

  # Loop through regions and perform analysis
  results_list <- list()
  
  for (i in seq_along(regions)) {
    region <- regions[i]
    
    tryCatch({
      log_line(sprintf("[%d/%d] Processing region: %s\n", i, length(regions), region))
      log_line(sprintf("[%s] Processing region: %s\n", Sys.time(), region))

      # Generate CTD for the region
      ctd <- generate_ctd(adata_prefix, region, cell_type_assignment_key)
      log_line(sprintf("[%s] CTD generated\n", Sys.time()))

      # Perform cell type enrichment analysis
      enrichment_results <- perform_cell_type_enrichment(gwas_base, ctd, region, ref_tax_output_dir, upstream_kb, downstream_kb, results_prefix)
      log_line(sprintf("[%s] Enrichment analysis complete\n", Sys.time()))

      # Plot results
      plot_results(enrichment_results, ctd, region, ref_tax_output_dir, gwas_base)
      log_line(sprintf("[%s] Plotting complete\n", Sys.time()))

      # Save results
      saveRDS(enrichment_results, file.path(ref_tax_output_dir, paste0("enrichment_results_", region, ".rds")))
      write.table(enrichment_results, file.path(ref_tax_output_dir, paste0("enrichment_results_", region, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
      log_line(sprintf("[%s] Results saved\n", Sys.time()))
      
      # Clean up
      rm(ctd, enrichment_results)
      gc()
      
      results_list[[region]] <- "SUCCESS"
      
    }, error = function(e) {
      log_line(sprintf("[%s] ERROR: %s\n", Sys.time(), e$message))
      log_line(sprintf("ERROR processing region %s: %s\n", region, e$message))
      results_list[[region]] <<- paste("FAILED:", e$message)
    })
  }

  log_line("Analysis complete\n")
  log_line(sprintf("Successful regions: %d\n", sum(results_list == "SUCCESS")))
  log_line(sprintf("Failed regions: %d\n", sum(results_list != "SUCCESS")))
}

if (!interactive()) {

  parser <- ArgumentParser(description = "Cell Type Enrichment Analysis")
  parser$add_argument("--regions", type = "character", required = TRUE, help = "space-separated list of regions to analyze")
  parser$add_argument("--data-dir", type = "character", default = NULL, help = "Data directory")
  parser$add_argument("--ref-tax-output-dir", type = "character", required = TRUE, help = "Reference taxonomy output directory")
  parser$add_argument("--gwas-base", type = "character", required = TRUE, help = "GWAS summary statistics file")
  parser$add_argument("--magma-gwas-base", type = "character", default = NULL, help = "Base name for MAGMA-formatted GWAS summary statistics file")
  parser$add_argument("--adata-prefix", type = "character", required = TRUE, help = "Prefix for anndata files")
  parser$add_argument("--workflow-mode", type = "character", required = TRUE, help = "Workflow mode: precomputed, full, or aibs_test")
  parser$add_argument("--upstream-kb", type = "integer", default = 35, help = "Upstream kb for SNP to gene mapping")
  parser$add_argument("--downstream-kb", type = "integer", default = 10, help = "Downstream kb for SNP to gene mapping")
  parser$add_argument("--cell-type-assignment-key", type = "character", required = TRUE, help = "Key in anndata.obs for cell type assignments")
  parser$add_argument("--results-prefix", type = "character", required = TRUE, help = "Prefix for results files")

  args <- parser$parse_args()

  regions <- trimws(unlist(strsplit(args$regions, " ")))
  workflow_mode <- args$workflow_mode
  ref_tax_output_dir <- args$ref_tax_output_dir
  adata_prefix <- args$adata_prefix
  gwas_base <- args$gwas_base
  upstream_kb <- args$upstream_kb
  downstream_kb <- args$downstream_kb
  cell_type_assignment_key <- args$cell_type_assignment_key
  results_prefix <- args$results_prefix


  if (workflow_mode == "full") {
    magma_gwas_base <- args$magma_gwas_base
    data_dir <- args$data_dir
    gwas_summary_stats <- file.path(data_dir, gwas_base, ".tsv")
    magma_gwas_summary_stats <- file.path(ref_tax_output_dir, magma_gwas_base, ".tsv")
  }

  if (workflow_mode == "precomputed") {
    main_precomputed(ref_tax_output_dir, gwas_base, adata_prefix, regions, upstream_kb, downstream_kb, cell_type_assignment_key, results_prefix)
  } else if (workflow_mode == "full") {
    main_full(ref_tax_output_dir, gwas_summary_stats, magma_gwas_summary_stats, adata_prefix, regions, upstream_kb, downstream_kb, cell_type_assignment_key, results_prefix)
  } else {
    stop("Invalid workflow mode")
  }
}