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
REDUCED_ADATA_PREFIX = config["reduced_adata_prefix"]
BATCH_KEY = config["batch_key"]
LATENT_KEY_SCVI = config["latent_key_scvi"]
MMC_CELL_TYPE_LABEL_KEY = config["mmc_cell_type_label_key"]
SCVI_ADATA_PREFIX = config["scvi_adata_prefix"]
LATENT_KEY_SCANVI = config["latent_key_scanvi"]
OUTPUT_SCVI_DIR = config["output_scvi_dir"]
SCANVI_ADATA_PREFIX = config["scanvi_adata_prefix"]
PREDICTIONS_KEY = config["scanvi_predictions_key"]
OUTPUT_SCANVI_DIR = config["output_scanvi_dir"]
OUTPUT_CELL_TYPES_FILE_PREFIX = config["scanvi_output_cell_types_file_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{SCVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

rule integrate_scvi:
    """Integrate data for each brain region using scVI"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{REDUCED_ADATA_PREFIX}",
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{SCVI_ADATA_PREFIX}",
        batch_key=BATCH_KEY,
        latent_key=LATENT_KEY_SCVI,
        output_scvi_dir=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_SCVI_DIR}"

    threads:
        4

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        gpu=4,
        gpu_type="nvidia-tesla-t4",
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/integrate_scvi.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run the integration script
        python3 -u GP2-Expansion/workflows/src/integrate_scvi.py \
            --adata-input-prefix {params.adata_input_prefix} \
            --adata-output-prefix {params.adata_output_prefix} \
            --output-scvi-dir {params.output_scvi_dir} \
            --batch-key {params.batch_key} \
            --latent-key {params.latent_key} \
            --regions '{params.regions}' \
            2>&1 | tee {log}
        """

rule label_scanvi:
    """Leverage cell-type from MMC to assign the rest of the cells with scANVI"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{SCVI_ADATA_PREFIX}",
        scvi_outputs_dir_prefix=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_SCVI_DIR}",
        scanvi_outputs_dir_prefix=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_SCANVI_DIR}",
        latent_key=LATENT_KEY_SCANVI,
        cell_type_label_key=MMC_CELL_TYPE_LABEL_KEY,
        predictions_key=PREDICTIONS_KEY,
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}",
        output_cell_types_file_prefix=f"{REF_TAX_OUTPUT_DIR}/{OUTPUT_CELL_TYPES_FILE_PREFIX}"

    threads:
        4

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        gpu=4,
        gpu_type="nvidia-tesla-t4",
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/assign_remaining_cells.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base-gpu:latest"

    shell:
        """
        set -euo pipefail

        # Run the label prediction script
        python3 -u GP2-Expansion/workflows/src/label_scanvi.py \
            --adata-input-prefix {params.adata_input_prefix} \
            --scvi-outputs-dir-prefix {params.scvi_outputs_dir_prefix} \
            --scanvi-outputs-dir-prefix {params.scanvi_outputs_dir_prefix} \
            --latent-key-scanvi {params.latent_key} \
            --cell-type-label-key {params.cell_type_label_key} \
            --predictions-key {params.predictions_key} \
            --adata-output-prefix {params.adata_output_prefix} \
            --output-cell-types-file-prefix {params.output_cell_types_file_prefix} \
            --regions "{params.regions}" \
            2>&1 | tee {log}
        """