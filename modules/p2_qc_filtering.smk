# Imports
import os
import math
import sys
from pathlib import Path

# Get absolute path to workflows directory in order to import helper functions
workflows_root = "GP2-Expansion/workflows"
sys.path.insert(0, str(workflows_root))
from src.helpers import parse_regions_file, total_bytes

# Configuration 
configfile: f"{workflows_root}/config.yml"

STANDALONE = config["standalone"]
CONTAINER_REGISTRY = config["container_registry"]
OUTPUT_DIR = config["output_dir"]
SAMPLE_ADATA_FILES = config["sample_adata_files"]
MERGED_ADATA_PREFIX = config["merged_adata_prefix"]
QC_PLOTS_DIR = config["qc_plots_dir"]
REGIONS_FILE = config["regions_file"]
QC_METRICS = config["qc_metrics"]
PCT_COUNTS_MT_MAX = config["pct_counts_mt_max"]
DOUBLETS_MAX = config["doublet_score_max"]
TOTAL_COUNTS_LIMITS = config["total_counts_limits"]
N_GENES_BY_COUNTS_LIMITS = config["n_genes_by_counts_limits"]
FILTERED_ADATA_PREFIX = config["filtered_adata_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{OUTPUT_DIR}/{QC_PLOTS_DIR}/violin_{{metric}}_{{region}}.png", metric=QC_METRICS.split(", "), region=REGIONS),
            expand(f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

rule plot_qc_metrics:
    """Merge preprocessed samples and generate QC plots"""
    input:
        expand(f"{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{OUTPUT_DIR}/{QC_PLOTS_DIR}/violin_{{metric}}_{{region}}.png", metric=QC_METRICS.split(", "), region=REGIONS),

    params:
        regions = " ".join(REGIONS),
        adata_input_prefix=f"{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}",
        qc_metrics=QC_METRICS,
        qc_plots_dir=f"{OUTPUT_DIR}/{QC_PLOTS_DIR}"

    threads: 
        12

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{OUTPUT_DIR}/logs/plot_qc_metrics.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail
        
        # Run plot QC
        python3 -u GP2-Expansion/workflows/src/plot_qc.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --qc-plots-dir {params.qc_plots_dir} \
            --metrics "{params.qc_metrics}" \
            2>&1 | tee {log}
        """

rule filter_adata:
    """Filter merged AnnData object based on QC metrics and doublet scores"""
    input:
        expand(f"{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    params:
        regions = " ".join(REGIONS),
        adata_input_prefix=f"{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}",
        pct_counts_mt_max=PCT_COUNTS_MT_MAX,
        doublet_score_max=DOUBLETS_MAX,
        total_counts_limits=TOTAL_COUNTS_LIMITS,
        n_genes_by_counts_limits=N_GENES_BY_COUNTS_LIMITS,
        adata_output_prefix=f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}"

    threads:
        1

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{OUTPUT_DIR}/logs/filter_adata.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run filtering
        python3 -u GP2-Expansion/workflows/src/filter.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --pct-counts-mt-max {params.pct_counts_mt_max} \
            --doublet-score-max {params.doublet_score_max} \
            --total-counts-limits {params.total_counts_limits} \
            --n-genes-by-counts-limits {params.n_genes_by_counts_limits} \
            --adata-output-prefix {params.adata_output_prefix} \
            2>&1 | tee {log}
        """
