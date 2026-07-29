# Imports
import os
from pathlib import Path
import pandas as pd
import yaml

# Helper function to load config from YAML file
def load_config() -> dict:
    """Load configuration from YAML file."""
    workflows_root = 'GP2-Expansion/workflows'
    config_path = f"{workflows_root}/config.yml"
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

# Helper functions to parse regions file and output space-separated regions list
def parse_regions_file(file_path: str) -> list:
    """
    Parse a text file containing regions (one per line) into a list.
    
    Args:
        file_path: Path to the text file containing regions.
    Returns:
        List of regions as strings.
    """
    regions = []
    try:
        with open(file_path, 'r') as f:
            regions = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading regions file {file_path}: {e}")
    return regions

# Helper functions to parse regions file and output space-separated regions list
def parse_datasets_file(file_path: str) -> list:
    """
    Parse a text file containing datasets (one per line) into a list.
    
    Args:
        file_path: Path to the text file containing datasets.
    Returns:
        List of datasets as strings.
    """
    datasets = []
    try:
        with open(file_path, 'r') as f:
            datasets = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading datasets file {file_path}: {e}")
    return datasets

# Helper function to normalize region names to eliminate inconsistencies in metadata
def normalize_region_names(region_series: pd.Series) -> pd.Series:
    """
    Normalize region names by stripping whitespace and converting to title case.
    
    Args:
        region_series: pandas Series with region names

    Returns:
        pandas Series with normalized region names.
    """
    return region_series.str.strip().str.title().str.replace(' ', '_')

# Helper function to calculate resources based on input file size
def total_bytes_from_file(file_path):
    """Get size of a single file in bytes."""
    if not os.path.exists(file_path):
        return 0
    try:
        return Path(file_path).stat().st_size
    except Exception as e:
        print(f"Warning: Could not get file size for {file_path}: {e}")
        return 0

# Helper function to calculate total size in bytes for a list of files
def total_bytes(
        adata_prefix = None,
        suffix = None,
        additional_files=None
        ):
    """
    Get size in bytes for local files.
    Args:
        adata_prefix: Prefix for AnnData files (e.g. "output/region_processed")
        suffix: List of suffixes to check for corresponding AnnData files
        additional_files: List of additional files to include in size calculation
    Returns:
        Total size in bytes of the specified files.
    """
    total_size = 0
    if not suffix is None and not adata_prefix is None:
        for s in suffix:
            adata_object = f"{adata_prefix}{s}.h5ad"
            file_size = total_bytes_from_file(adata_object)
            total_size += file_size
    
    if not additional_files is None:
        for file in additional_files:
            file_size = total_bytes_from_file(file)
            total_size += file_size

    return total_size

def total_bytes_largest_region(
        sample_adata_files: str
        ) -> list:
    """
    Get list of file paths for the largest sample AnnData files for a specific region.
    Args:
        sample_adata_files: Path to the file listing sample adata files (tab-separated).
    Returns:
        List of file paths for the largest sample AnnData files for the specified region.
    """
    sample_adata_files_df = pd.read_csv(sample_adata_files, sep='\t')
    config = load_config()
    region_col = config['metadata_region_col']
    adata_file_col = config['adata_file_col']
    regions = sample_adata_files_df[region_col].unique()
    largest_size = 0
    for region in regions:
        region_sample_files_df = sample_adata_files_df[sample_adata_files_df[region_col] == region]
        region_size = 0
        region_file_paths = region_sample_files_df[adata_file_col].tolist()
        region_size = total_bytes(additional_files=region_file_paths)
        if region_size > largest_size:
            largest_size = region_size
    
    return largest_size
