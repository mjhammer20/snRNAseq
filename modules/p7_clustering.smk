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
LATENT_KEY = config["latent_key_scvi"]
SCANVI_ADATA_PREFIX = config["scanvi_adata_prefix"]
N_NEIGHBORS = config["clustering_n_neighbors"]
LEIDEN_RES = config["leiden_res"]
CLUSTERED_ADATA_PREFIX = config["clustered_adata_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

rule clustering:
    """Cluster the latent space with UMAP and Leiden clustering"""
    input:
        expand(f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{REF_TAX_OUTPUT_DIR}/{SCANVI_ADATA_PREFIX}",
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{CLUSTERED_ADATA_PREFIX}",
        latent_key=LATENT_KEY,
        n_neighbors=N_NEIGHBORS,
        leiden_res=LEIDEN_RES

    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/clustering.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run clustering script
        python3 -u GP2-Expansion/workflows/src/clustering_umap.py \
            --adata-input-prefix {params.adata_input_prefix} \
            --adata-output-prefix {params.adata_output_prefix} \
            --latent-key {params.latent_key} \
            --n-neighbors {params.n_neighbors} \
            --leiden-res {params.leiden_res} \
            --regions "{params.regions}" \
            2> {log}
        """


