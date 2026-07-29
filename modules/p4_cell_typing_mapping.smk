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
REF_TAX = config["mmc_ref_tax"]
REF_TAX_OUTPUT_DIR = f"{OUTPUT_DIR}/{REF_TAX}"
REF_DIR = config["mmc_ref_dir"]
REF_PATH = f'{DATA_DIR}/{REF_DIR}'
MMC_GENE_MAPPER_DIR = config["mmc_gene_mapper_dir"]
GENE_MAPPER_PATH = f'{DATA_DIR}/{MMC_GENE_MAPPER_DIR}'
MMC_MARKER_GENES = config["mmc_marker_genes"]
MMC_PRECOMPUTED_STATS = config["mmc_precomputed_stats"]
MMC_GENE_MAPPER_DB = config["mmc_gene_mapper_db"]
META_ANNOTATED_ADATA_PREFIX = config["meta_annotated_adata_prefix"]
MMC_RESULTS_PREFIX = config["mmc_results_prefix"]
N_PROCESSORS = config["mmc_n_processors"]
MAX_GB = config["mmc_max_gb"]
RNG_SEED = config["mmc_rng_seed"]
MMC_CHUNK_SIZE = config["mmc_chunk_size"]
N_RUNNERS_UP = config["mmc_n_runners_up"]
MMC_ANNOTATED_ADATA_PREFIX = config["mmc_annotated_adata_prefix"]
MMC_OUTPUT_CELL_TYPES_FILE_PREFIX = config["mmc_output_cell_types_file_prefix"]

# Parse regions file to get list of regions for resource calculation and output file naming
REGIONS = parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")

# Compute total GB for merged adata files once at parse time for resource calculation
total_gb = total_bytes(adata_prefix=f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}_", suffix=REGIONS) / (1024 ** 3)

# Create reference taxonomy output directory, and associated logs directory, if it doesn't exist
Path(REF_TAX_OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(f"{REF_TAX_OUTPUT_DIR}/logs").mkdir(parents=True, exist_ok=True)

# Rules

if STANDALONE:
    rule all:
        input:
            expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.extended_results.json", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.results.csv", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.log.txt", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
            expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_OUTPUT_CELL_TYPES_FILE_PREFIX}_{{region}}.parquet", region=REGIONS)

rule run_mmc:
    """Run Cell Type Mapper for cell type annotation"""
    input:
        mmc_marker_genes=f'{REF_PATH}/{MMC_MARKER_GENES}',
        mmc_precomputed_stats=f'{REF_PATH}/{MMC_PRECOMPUTED_STATS}',
        mmc_gene_mapper_db=f'{GENE_MAPPER_PATH}/{MMC_GENE_MAPPER_DB}'

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.extended_results.json", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.results.csv", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.log.txt", region=REGIONS)

    params:
        n_processors=N_PROCESSORS,
        max_gb=MAX_GB,
        rng_seed=RNG_SEED,
        chunk_size=MMC_CHUNK_SIZE,
        n_runners_up=N_RUNNERS_UP,
        adata_prefix=f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}",
        output_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}",
        regions=" ".join(REGIONS)

    threads:
        8

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/run_mmc.log"

    container:
        f"{CONTAINER_REGISTRY}/gp2-omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run MMC
        python3 -u GP2-Expansion/workflows/src/mmc.py \
            --adata-input-prefix {params.adata_prefix} \
            --regions "{params.regions}" \
            --mmc-marker-genes {input.mmc_marker_genes} \
            --mmc-precomputed-stats {input.mmc_precomputed_stats} \
            --mmc-gene-mapper-db {input.mmc_gene_mapper_db} \
            --n-processors {params.n_processors} \
            --max-gb {params.max_gb} \
            --rng-seed {params.rng_seed} \
            --chunk-size {params.chunk_size} \
            --n-runners-up {params.n_runners_up} \
            --output-prefix {params.output_prefix} \
            2>&1 | tee {log}
        """
    
rule add_cell_type_annotations:
    """Annotate cell types for each brain region using transcriptional phenotype workflow"""
    input:
        expand(f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}_{{region}}.results.csv", region=REGIONS)

    output:
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}_{{region}}.h5ad", region=REGIONS),
        expand(f"{REF_TAX_OUTPUT_DIR}/{MMC_OUTPUT_CELL_TYPES_FILE_PREFIX}_{{region}}.parquet", region=REGIONS)

    params:
        output_dir=OUTPUT_DIR,
        regions=" ".join(REGIONS),
        adata_input_prefix=f"{OUTPUT_DIR}/{META_ANNOTATED_ADATA_PREFIX}",
        mmc_results_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_RESULTS_PREFIX}",
        adata_output_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_ANNOTATED_ADATA_PREFIX}",
        output_cell_types_file_prefix=f"{REF_TAX_OUTPUT_DIR}/{MMC_OUTPUT_CELL_TYPES_FILE_PREFIX}"

    threads:
        4

    resources:
        mem_mb = lambda wildcards, input: max(int((total_gb * 18 + 20) * 1024), 131072),
        disk_mb = lambda wildcards, input: max(int((total_gb * 2 + 5) * 1024), 10240),
        runtime=10800

    log:
        f"{REF_TAX_OUTPUT_DIR}/logs/annotate_cell_types.log"

    container:
        f"{CONTAINER_REGISTRY}/omics-base:latest"

    shell:
        """
        set -euo pipefail

        # Run transcriptional phenotype annotation script
        python3 -u GP2-Expansion/workflows/src/transcriptional_phenotype.py \
            --adata-input-prefix {params.adata_input_prefix} \
            --mmc-results-prefix {params.mmc_results_prefix} \
            --regions "{params.regions}" \
            --adata-output-prefix {params.adata_output_prefix} \
            --output-cell-types-file-prefix {params.output_cell_types_file_prefix} \
            2>&1 | tee {log}
        """