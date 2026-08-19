# Imports
import os
import math
import sys
from pathlib import Path

# Get absolute path to workflows directory in order to import helper functions
WORKFLOW_ROOT = str(
    _snakefile_path.parent.parent
    if _snakefile_path.parent.name == "modules"
    else _snakefile_path.parent
)
sys.path.insert(0, str(WORKFLOW_ROOT))
from src.helpers import parse_regions_file, total_bytes

# Configuration 
configfile: f"{WORKFLOW_ROOT}/config.yml"

CONTAINER = config["container"]
STANDALONE = config["standalone"]
OUTPUT_DIR = config["output_dir"]
REGIONS_FILE = config["regions_file"]
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
CLUSTERED_ADATA_PREFIX = config["clustered_adata_prefix"]
TARGET_SUM = config["target_sum"]
MMC_CELL_TYPE_LABEL_KEY = config["mmc_cell_type_label_key"]
SCANVI_PREDICTION_KEY = config["scanvi_predictions_key"]
CELL_TYPE_ASSIGNMENT_KEY = config["cell_type_assignment_key"]
FINAL_ADATA_PREFIX = config["final_adata_prefix"]
FINAL_ADATA_METADATA_PREFIX = config["final_adata_metadata_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
TOTAL_GB = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_METADATA_PREFIX}_{{region}}.tsv", region=REGIONS)

rule prepare_final_adata:
    """Assign cell types based on cell type mappings from MMC (priority) and scANVI predictions (backup)"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_METADATA_PREFIX}_{{region}}.tsv", region=REGIONS)

    params:
        workflow_root=WORKFLOW_ROOT,
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}",
        target_sum=TARGET_SUM,
        mmc_cell_type_label_key=MMC_CELL_TYPE_LABEL_KEY,
        scanvi_prediction_key=SCANVI_PREDICTION_KEY,
        cell_type_assignment_key=CELL_TYPE_ASSIGNMENT_KEY,
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_PREFIX}",
        metadata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{FINAL_ADATA_METADATA_PREFIX}"

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((TOTAL_GB * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((TOTAL_GB * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/cell_type_assignment.log"

    container:
        CONTAINER

    shell:
        """
        set -euo pipefail

        # Run cell type assignment script
        python3 -u {params.workflow_root}/src/prepare_final_adata.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --target-sum {params.target_sum} \
            --mmc-cell-type-label-key "{params.mmc_cell_type_label_key}" \
            --scanvi-prediction-key "{params.scanvi_prediction_key}" \
            --cell-type-assignment-key "{params.cell_type_assignment_key}" \
            --adata-output-prefix {params.adata_output_prefix} \
            --metadata-output-prefix {params.metadata_output_prefix} \
            2>&1 | tee {log}
        """

