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

CONTAINER_REGISTRY = config["container_registry"]
STANDALONE = config["standalone"]
OUTPUT_DIR = config["output_dir"]
REGIONS_FILE = config["regions_file"]
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
FINAL_ADATA_PREFIX = config["final_adata_prefix"]
DESEQ2_RESULTS_PREFIX = config["deseq2_results_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{DESEQ2_RESULTS_PREFIX}_{{region}}.tsv", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{DESEQ2_RESULTS_PREFIX}_all_significant_genes.tsv")


rule identify_differentially_expressed_genes:
    """Identify differentially expressed genes for each region"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{DESEQ2_RESULTS_PREFIX}_{{region}}.tsv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{DESEQ2_RESULTS_PREFIX}_all_significant_genes.tsv")

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        results_prefix=f"{REF_TAX_OUTPUT_DIR}/{DESEQ2_RESULTS_PREFIX}"

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/differential_expression_analysis.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail
        
        # Run differential expression script
        python3 -u GP2-Expansion/workflows/src/differential_expression.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --results-prefix {params.results_prefix} \
            2>&1 | tee {log}
        """