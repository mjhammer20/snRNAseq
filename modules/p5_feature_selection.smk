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
OUTPUT_DIR = config["output_dir"]
DATA_DIR = config["data_dir"]
REGIONS_FILE = config["regions_file"]
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
ENSEMBL_MAPPING_FILE = config["ensembl_mapping_file"]
GWAS_SUMMARY_STATS = config["gwas_summary_stats"]
MMC_RESULTS_PREFIX = config["mmc_results_prefix"]
MMC_ANNOTATED_ADATA_PREFIX = config["mmc_annotated_adata_prefix"]
BATCH_KEY = config["batch_key"]
N_TOP_GENES = config["n_top_genes"]
REDUCED_ADATA_PREFIX = config["reduced_adata_prefix"]
OUTPUT_ALL_GENES_PREFIX = config["output_all_genes_prefix"]
OUTPUT_HVG_GENES_PREFIX = config["output_hvg_genes_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
rule all:
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_ALL_GENES_PREFIX}_{{region}}.csv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_HVG_GENES_PREFIX}_{{region}}.csv", region=REGIONS)

rule process:
    """Process filtered AnnData object for each brain region"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.extended_results.json", region=REGIONS),
        gwas_summary_stats=f"{DATA_DIR}/{GWAS_SUMMARY_STATS}",
        gene_mapping_file=f"{DATA_DIR}/{ENSEMBL_MAPPING_FILE}"

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_ALL_GENES_PREFIX}_{{region}}.csv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_HVG_GENES_PREFIX}_{{region}}.csv", region=REGIONS)

    params:
        mmc_results_prefix=f'{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}',
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}",
        batch_key=BATCH_KEY,
        n_top_genes=N_TOP_GENES,
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}",
        output_all_genes_prefix=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_ALL_GENES_PREFIX}",
        output_hvg_genes_prefix=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_HVG_GENES_PREFIX}", 
        regions=" ".join(REGIONS)

    threads:
        4

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/feature_selection.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        python3 -u GP2-Expansion/workflows/src/feature_selection.py \
            --adata-input-prefix {params.adata_input_prefix} \
            --regions "{params.regions}" \
            --gwas-summary-stats {input.gwas_summary_stats} \
            --mmc-results-prefix {params.mmc_results_prefix} \
            --gene-mapping-file {input.gene_mapping_file} \
            --batch-key {params.batch_key} \
            --n-top-genes {params.n_top_genes} \
            --adata-output-prefix {params.adata_output_prefix} \
            --output-all-genes-prefix {params.output_all_genes_prefix} \
            --output-hvg-genes-prefix {params.output_hvg_genes_prefix} \
            2>&1 | tee {log}
        """