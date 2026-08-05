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
DATA_DIR = config["data_dir"]
REGIONS_FILE = config["regions_file"]
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
FINAL_ADATA_PREFIX = config["final_adata_prefix"]
GWAS_SUMMARY_STATS = config["gwas_summary_stats"]
CELL_TYPE_ASSIGNMENT_KEY = config["cell_type_assignment_key"]
CONDITION_KEY = config["condition_key"]
GENE_SYMBOLS_KEY = config["gene_symbols_key"]
MIN_CELLS = config["min_cells"]
MEAN_EXPRESSION_THRESHOLD = config["mean_expression_threshold"]
MEAN_FRACTION_EXPRESSED_THRESHOLD = config["mean_fraction_expressed_threshold"]
RESULTS = config["results"]
RESULTS_DIR = f"{REF_TAX_OUTPUT_DIR}/{RESULTS}"
DOTPLOTS_PREFIX = config["dotplots_prefix"]
EXPRESSED_GENES_PREFIX = config["expressed_genes_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
rule all:
    input:
        expand(f"{RESULTS_DIR}/{DOTPLOTS_PREFIX}_{{region}}.png", region=REGIONS),
        expand(f"{RESULTS_DIR}/{EXPRESSED_GENES_PREFIX}_{{region}}.tsv", region=REGIONS)

rule cell_type_expression_analysis:
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        gwas_summary_stats_path = f"{DATA_DIR}/{GWAS_SUMMARY_STATS}"

        
    output:
        expand(f"{RESULTS_DIR}/{DOTPLOTS_PREFIX}_{{region}}.png", region=REGIONS),
        expand(f"{RESULTS_DIR}/{EXPRESSED_GENES_PREFIX}_{{region}}.tsv", region=REGIONS)

    params:
        regions=REGIONS,
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        cell_type_assignment_key=CELL_TYPE_ASSIGNMENT_KEY,
        gene_symbols_key=GENE_SYMBOLS_KEY,
        condition_key=CONDITION_KEY,
        min_cells=MIN_CELLS,
        mean_expression_threshold=MEAN_EXPRESSION_THRESHOLD,
        mean_fraction_expressed_threshold=MEAN_FRACTION_EXPRESSED_THRESHOLD,
        results_dir=RESULTS_DIR,
        expressed_genes_prefix=EXPRESSED_GENES_PREFIX,
        dotplots_prefix=DOTPLOTS_PREFIX

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/cell_type_expression_analysis.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        python3 -u GP2-Expansion/workflows/src/cell_type_expression_analysis.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --gwas-summary-stats {input.gwas_summary_stats_path} \
            --cell-type-assignment-key "{params.cell_type_assignment_key}" \
            --gene-symbols-key "{params.gene_symbols_key}" \
            --condition-key "{params.condition_key}" \
            --min-cells {params.min_cells} \
            --mean-threshold {params.mean_expression_threshold} \
            --frac-threshold {params.mean_fraction_expressed_threshold} \
            --results-dir {params.results_dir} \
            --expressed-genes-prefix {params.expressed_genes_prefix} \
            --dotplots-prefix {params.dotplots_prefix} \
            2>&1 | tee {log}
        """