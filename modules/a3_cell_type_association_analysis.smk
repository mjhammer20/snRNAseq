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
CELL_TYPE_ASSIGNMENT_KEY = config["cell_type_assignment_key"]
WORKFLOW_MODE = config["workflow_mode"]
GWAS_BASE = config["gwas_base"]
MAGMA_GWAS_BASE = config["magma_gwas_base"]
UPSTREAM_KB = config["upstream_kb"]
DOWNSTREAM_KB = config["downstream_kb"]
MAGMA_RESULTS_PREFIX = config["magma_results_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Determine analysis name based on workflow mode
if WORKFLOW_MODE == "precomputed":
    ANALYSIS_NAME = GWAS_BASE
else:
    ANALYSIS_NAME = MAGMA_GWAS_BASE

# Rules
rule all:
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{MAGMA_RESULTS_PREFIX}_{{region}}.tsv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/MAGMA_Figures/{ANALYSIS_NAME}/{ANALYSIS_NAME}.{UPSTREAM_KB}UP.{DOWNSTREAM_KB}DOWN.annotLevel1.ConditionalFacets.{{region}}.Merged.pdf", region=REGIONS)


rule cell_type_association_analysis:
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)
        
    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{MAGMA_RESULTS_PREFIX}_{{region}}.tsv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/MAGMA_Figures/{ANALYSIS_NAME}/{ANALYSIS_NAME}.{UPSTREAM_KB}UP.{DOWNSTREAM_KB}DOWN.annotLevel1.ConditionalFacets.{{region}}.Merged.pdf", region=REGIONS)

    params:
        regions=REGIONS,
        data_dir=DATA_DIR,
        ref_tax_output_dir=REF_TAX_OUTPUT_DIR,
        gwas_base=GWAS_BASE,
        magma_gwas_base=MAGMA_GWAS_BASE,
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        workflow_mode=WORKFLOW_MODE,
        upstream_kb=UPSTREAM_KB,
        downstream_kb=DOWNSTREAM_KB,
        cell_type_assignment_key=CELL_TYPE_ASSIGNMENT_KEY,
        results_prefix=MAGMA_RESULTS_PREFIX

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/cell_type_association_analysis.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base-R:latest"

    shell:
        """
        Rscript GP2-Expansion/workflows/src/cell_type_enrichment.r \
            --regions "{params.regions}" \
            --data-dir {params.data_dir} \
            --ref-tax-output-dir {params.ref_tax_output_dir} \
            --gwas-base {params.gwas_base} \
            --magma-gwas-base {params.magma_gwas_base} \
            --adata-prefix {params.adata_input_prefix} \
            --workflow-mode {params.workflow_mode} \
            --upstream-kb {params.upstream_kb} \
            --downstream-kb {params.downstream_kb} \
            --cell-type-assignment-key "{params.cell_type_assignment_key}" \
            --results-prefix {params.results_prefix} \
            2>&1 | tee {log}
        """