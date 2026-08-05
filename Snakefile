# ============================================================
# snRNAseq Full Workflow Snakefile
# Runs all processing (p1–p9) and analysis (a1–a3) modules
# in order. Set standalone: FALSE in config.yml when using
# this master Snakefile (module-level `rule all` blocks are
# guarded by the STANDALONE flag and will not conflict).
# ============================================================

import os
import sys
from pathlib import Path

# Resolve the workflow root so module-level helper imports work
workflows_root = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.insert(0, workflows_root)

from src.helpers import parse_regions_file

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
configfile: "config.yml"

OUTPUT_DIR          = config["output_dir"]
REGIONS_FILE        = config["regions_file"]
REF_TAX             = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR  = f"{OUTPUT_DIR}/{REF_TAX}"

# REGIONS is produced at runtime by p1. All downstream expand()
# calls use a lambda so they resolve lazily after p1 completes.
def regions(wildcards):
    regions_path = checkpoints.list_sample_adata_files.get(**wildcards).output.regions
    return parse_regions_file(str(regions_path))

# ------------------------------------------------------------
# Include all modules (standalone guards prevent duplicate
# `rule all` definitions when standalone: FALSE)
# ------------------------------------------------------------

# Processing
include: "modules/p1_sample_merging.smk"
include: "modules/p2_qc_filtering.smk"
include: "modules/p3_sample_metadata_annotation.smk"
include: "modules/p4_cell_typing_mapping.smk"
include: "modules/p5_feature_selection.smk"
include: "modules/p6_integration.smk"
include: "modules/p7_clustering.smk"
include: "modules/p8_final_adata_prep.smk"
include: "modules/p9_report_generation.smk"

# Analysis
include: "modules/a1_cell_type_expression_analysis.smk"
include: "modules/a2_differential_expression_analysis.smk"
include: "modules/a3_cell_type_association_analysis.smk"

# ------------------------------------------------------------
# Master rule — collects all terminal outputs from every module
# ------------------------------------------------------------
rule all:
    input:
        # ── p1: Sample Merging ────────────────────────────────
        f"{OUTPUT_DIR}/{config['sample_adata_files']}",
        f"{OUTPUT_DIR}/{config['regions_file']}",
        lambda wc: expand(
            f"{OUTPUT_DIR}/{config['merged_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{OUTPUT_DIR}/{config['initial_adata_metadata_prefix']}_{{region}}.csv",
            region=regions(wc)
        ),

        # ── p2: QC Filtering ─────────────────────────────────
        lambda wc: expand(
            f"{OUTPUT_DIR}/{config['qc_plots_dir']}/violin_{{metric}}_{{region}}.png",
            metric=config["qc_metrics"].split(", "),
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{OUTPUT_DIR}/{config['filtered_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),

        # ── p3: Sample Metadata Annotation ───────────────────
        lambda wc: expand(
            f"{OUTPUT_DIR}/{config['meta_annotated_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),

        # ── p4: Cell Typing Mapping ───────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['mmc_results_prefix']}_{{region}}.extended_results.json",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['mmc_results_prefix']}_{{region}}.results.csv",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['mmc_annotated_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['mmc_output_cell_types_file_prefix']}_{{region}}.parquet",
            region=regions(wc)
        ),

        # ── p5: Feature Selection ─────────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['reduced_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['output_all_genes_prefix']}_{{region}}.csv",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['output_hvg_genes_prefix']}_{{region}}.csv",
            region=regions(wc)
        ),

        # ── p6: Integration ───────────────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['scvi_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['scanvi_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),

        # ── p7: Clustering ────────────────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['clustered_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),

        # ── p8: Final AnnData Prep ────────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['final_adata_prefix']}_{{region}}.h5ad",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['final_adata_metadata_prefix']}_{{region}}.tsv",
            region=regions(wc)
        ),

        # ── p9: Report Generation ─────────────────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['umap_dir']}/umap_features_{{region}}.png",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['umap_dir']}/umap_groups_{{region}}.png",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['scib_report_dir']}/scib_report_{{region}}.csv",
            region=regions(wc)
        ),

        # ── a1: Cell Type Expression Analysis ────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['results']}/{config['dotplots_prefix']}_{{region}}.png",
            region=regions(wc)
        ),
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['results']}/{config['expressed_genes_prefix']}_{{region}}.tsv",
            region=regions(wc)
        ),

        # ── a2: Differential Expression Analysis ─────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['deseq2_results_prefix']}_{{region}}.tsv",
            region=regions(wc)
        ),
        f"{REF_TAX_OUTPUT_DIR}/{config['deseq2_results_prefix']}_all_significant_genes.tsv",

        # ── a3: Cell Type Association Analysis ────────────────
        lambda wc: expand(
            f"{REF_TAX_OUTPUT_DIR}/{config['magma_results_prefix']}_{{region}}.tsv",
            region=regions(wc)
        ),
