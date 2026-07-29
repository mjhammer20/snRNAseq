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
DATA_DIR = config["data_dir"]
REGIONS_FILE = config["regions_file"]
SAMPLE_METADATA = config["sample_metadata"]
FILTERED_ADATA_PREFIX = config["filtered_adata_prefix"]
METADATA_SAMPLE_COL = config["metadata_sample_col"]
METADATA_REGION_COL = config["metadata_region_col"]
ADATA_SAMPLE_COL = config["adata_sample_col"]
META_ANNOTATED_ADATA_PREFIX = config["meta_annotated_adata_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Rules
if STANDALONE:
    rule all:
        input:
            expand(f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

rule add_sample_metadata:
    """Add sample metadata to region-specific AnnData objects"""
    input:
        expand(f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        sample_metadata=f"{DATA_DIR}/asap_sn_rna_metadata/{SAMPLE_METADATA}"

    output:
        expand(f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS)

    params:
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{OUTPUT_DIR}/{FILTERED_ADATA_PREFIX}",
        metadata_region_col=METADATA_REGION_COL,
        metadata_sample_col=METADATA_SAMPLE_COL,
        adata_sample_col=ADATA_SAMPLE_COL,
        adata_output_prefix=f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}"
        
    threads:
        16

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{OUTPUT_DIR}/logs/add_sample_metadata.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run sample metadata addition script
        python3 -u GP2-Expansion/workflows/src/add_sample_metadata.py \
            --regions "{params.regions}" \
            --adata-input-prefix {params.adata_input_prefix} \
            --metadata-region-col "{params.metadata_region_col}" \
            --sample-metadata {input.sample_metadata} \
            --metadata-sample-col "{params.metadata_sample_col}" \
            --adata-sample-col "{params.adata_sample_col}" \
            --adata-output-prefix {params.adata_output_prefix} \
            2>&1 | tee {log}
        """