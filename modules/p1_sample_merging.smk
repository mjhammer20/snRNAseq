# Imports
import os
import math
import sys
from pathlib import Path
import pandas as pd

# Get absolute path to workflows directory in order to import helper functions
workflows_root = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.insert(0, str(workflows_root))
from src.helpers import parse_regions_file, total_bytes_largest_region

# Configuration 
configfile: f"{workflows_root}/config.yml"

STANDALONE = config['standalone']
CONTAINER_REGISTRY = config['container_registry']
DATA_DIR = config['data_dir']
ADATA_DIR = config['adata_dir']
METADATA_DIR = config['metadata_dir']
OUTPUT_DIR = config['output_dir']
SAMPLE_METADATA = config['sample_metadata']
REGION_COL = config['metadata_region_col']
SAMPLE_COL = config['metadata_sample_col']
ADATA_FILE_COL = config['adata_file_col']
SAMPLE_ADATA_SUFFIX = config['sample_adata_suffix']
SAMPLE_ADATA_FILES = config['sample_adata_files']
REGIONS_FILE = config['regions_file']
MERGED_ADATA_PREFIX = config['merged_adata_prefix']
INITIAL_ADATA_METADATA_PREFIX = config['initial_adata_metadata_prefix']

# Rules
if STANDALONE:
    rule all:
        '''In standalone mode, the final output files are the merged AnnData objects and initial metadata files for each chunk.'''
        input:
            f'{OUTPUT_DIR}/{SAMPLE_ADATA_FILES}',
            f'{OUTPUT_DIR}/{REGIONS_FILE}',
            lambda wildcards: expand(f'{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}_{{region}}.h5ad', region=parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")),
            lambda wildcards: expand(f'{OUTPUT_DIR}/{INITIAL_ADATA_METADATA_PREFIX}_{{region}}.csv', region=parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}"))

rule list_sample_adata_files:
    '''List sample name, file path, and associated region in tsv file for downstream processing'''
    input:
        sample_metadata=f'{DATA_DIR}/{METADATA_DIR}/{SAMPLE_METADATA}'
    output:
        sample_adata_files=f'{OUTPUT_DIR}/{SAMPLE_ADATA_FILES}',
        regions = f'{OUTPUT_DIR}/{REGIONS_FILE}'
    params:
        adata_dir=f'{DATA_DIR}/{ADATA_DIR}',
        sample_adata_suffix=SAMPLE_ADATA_SUFFIX,
        sample_col=SAMPLE_COL,
        region_col=REGION_COL,
        adata_file_col=ADATA_FILE_COL
    log:
        f'{OUTPUT_DIR}/logs/list_sample_adata_files.log'
    threads: 1
    resources:
        mem_mb=256,
        disk_mb=10240
    container:
        f'{CONTAINER_REGISTRY}/omics-base:latest'
    shell:
        '''
        set -euo pipefail

        # Run list sample adata files script
        python3 -u GP2-Expansion/workflows/src/list_sample_adata_files.py \
            --sample-metadata {input.sample_metadata} \
            --adata-dir {params.adata_dir} \
            --sample-adata-suffix {params.sample_adata_suffix} \
            --sample-col "{params.sample_col}" \
            --region-col "{params.region_col}" \
            --adata-file-col "{params.adata_file_col}" \
            --output-sample-adata-files {output.sample_adata_files} \
            --output-regions {output.regions} \
            2>&1 | tee {log}
        '''


rule merge_sample_adatas:
    '''Merge preprocessed samples'''
    input:
        sample_adata_files=f'{OUTPUT_DIR}/{SAMPLE_ADATA_FILES}'
    output:
        expand(f'{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}_{{region}}.h5ad', region=parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}")),
        expand(f'{OUTPUT_DIR}/{INITIAL_ADATA_METADATA_PREFIX}_{{region}}.csv', region=parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}"))
    params:
        adata_output_prefix=f'{OUTPUT_DIR}/{MERGED_ADATA_PREFIX}',
        output_metadata_prefix=f'{OUTPUT_DIR}/{INITIAL_ADATA_METADATA_PREFIX}',
        region_col=REGION_COL,
        sample_col=SAMPLE_COL,
        adata_file_col=ADATA_FILE_COL,
        regions = " ".join(parse_regions_file(f"{OUTPUT_DIR}/{REGIONS_FILE}"))
    threads: 12
    resources:
        mem_mb=lambda wildcards, input: max(int(total_bytes_largest_region(input.sample_adata_files) * 10 / (1024 ** 2)), 131072), 
        disk_mb=lambda wildcards, input: max(int(total_bytes_largest_region(input.sample_adata_files) * 2 / (1024 ** 2)), 10240)
    log:
        f'{OUTPUT_DIR}/logs/merge_adata.log'
    container:
        f'{CONTAINER_REGISTRY}/omics-base:latest'
    shell:
        '''
        set -euo pipefail
        
        # Run merge adata script
        python3 -u GP2-Expansion/workflows/src/merge_samples_for_region.py \
            --regions "{params.regions}" \
            --sample-adata-files {input.sample_adata_files} \
            --sample-col "{params.sample_col}" \
            --region-col "{params.region_col}" \
            --adata-file-col "{params.adata_file_col}" \
            --adata-output-prefix {params.adata_output_prefix} \
            --output-metadata-prefix {params.output_metadata_prefix} \
            2>&1 | tee {log}
        '''