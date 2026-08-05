# Imports
import os
import math
import sys
from pathlib import Path

# Get absolute path to workflows directory in order to import helper functions
workflows_root = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.insert(0, str(workflows_root))
from src.helpers import parse_regions_file, total_bytes

# Configuration 
configfile: f"{workflows_root}/config.yml"

CONTAINER_REGISTRY = config["container_registry"]
STANDALONE = config["standalone"]
OUTPUT_DIR = config["output_dir"]
REGIONS_FILE = config["regions_file"]
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
FINAL_ADATA_PREFIX = config["final_adata_prefix"]
UMAP_GROUPS = config["umap_groups"]
UMAP_FEATURES = config["umap_features"]
UMAP_DIR = config["umap_dir"]
SCIB_REPORT_DIR = config["scib_report_dir"]
BATCH_KEY = config["batch_key"]
CELL_TYPE_ASSIGNMENT_KEY = config["cell_type_assignment_key"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{UMAP_DIR}/umap_features_{{region}}.png", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{UMAP_DIR}/umap_groups_{{region}}.png", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{SCIB_REPORT_DIR}/scib_report_{{region}}.csv", region=REGIONS)

rule generate_umaps:
    """Generate UMAP plots colored by groups and features for each region"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{UMAP_DIR}/umap_features_{{region}}.png", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{UMAP_DIR}/umap_groups_{{region}}.png", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        groups=UMAP_GROUPS,
        features=UMAP_FEATURES,
        umap_dir=f"{REF_TAX_OUTPUT_DIR}/{UMAP_DIR}"

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/generate_umaps.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run UMAP plotting script
        python3 -u GP2-Expansion/workflows/src/plot_groups_and_feats.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --umap-dir {params.umap_dir} \
            --groups "{params.groups}" \
            --features "{params.features}" \
            2>&1 | tee {log}
        """

rule calculate_artifact_metrics:
    """Calculate `scib` metrics for each region"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCIB_REPORT_DIR}/{region}/scib_report.csv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCIB_REPORT_DIR}/{region}/scib_results.svg", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        batch_key=BATCH_KEY,
        cell_type_assignment_key=CELL_TYPE_ASSIGNMENT_KEY,
        report_dir=f"{REF_TAX_OUTPUT_DIR}/{SCIB_REPORT_DIR}"

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/artifact_metrics.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run artifact metrics script
        python3 -u GP2-Expansion/workflows/src/artifact_metrics.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --batch-key {params.batch_key} \
            --cell-type-assignment-key {params.cell_type_assignment_key} \
            --output-report-dir {params.report_dir} \
            2>&1 | tee {log}
        """