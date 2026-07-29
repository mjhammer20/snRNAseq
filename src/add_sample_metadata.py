# Imports
import argparse
import scanpy as sc
from anndata import AnnData
import pandas as pd
import gc
import sys

# Get absolute path to workflows directory in order to import helper functions
workflows_root = 'GP2-Expansion/workflows'
sys.path.insert(0, str(workflows_root))
from src.helpers import normalize_region_names

def load_sample_metadata(file_path: str, metadata_region_col: str, metadata_sample_col: str) -> pd.DataFrame:
    """
    Load sample metadata from a comma-separated values (CSV) file.

    Args:
        file_path: Path to the CSV file containing sample metadata.
        metadata_region_col: Column name in sample_metadata that contains region information.
        metadata_sample_col: Column name in sample_metadata that contains sample identifiers.

    Returns:
        DataFrame containing sample metadata.
    """

    # Load sample metadata
    sample_metadata = pd.read_csv(file_path)

    # Only convert object columns to string dtype, don't fill yet
    for col in sample_metadata.columns:
        if sample_metadata[col].dtype == object:
            sample_metadata[col] = sample_metadata[col].astype(str)

    # Normalize region names in sample_metadata to ensure consistency
    sample_metadata[metadata_region_col] = normalize_region_names(sample_metadata[metadata_region_col])

    # Ensure sample column is string type for merging
    sample_metadata[metadata_sample_col] = sample_metadata[metadata_sample_col].astype(str)

    return sample_metadata

# Function to add sample metadata to the merged AnnData object
def add_sample_metadata(adata: AnnData, sample_metadata: pd.DataFrame, metadata_region_col: str, metadata_sample_col: str, adata_sample_col: str) -> AnnData:
    """
    Add sample metadata to the merged AnnData object.

    Args:
        adata: Merged AnnData object containing all samples.
        sample_metadata: DataFrame containing sample metadata.
        metadata_region_col: Column name in sample_metadata and adata.obs that contains region information.
        metadata_sample_col: Column name in sample_metadata that contains sample identifiers.
        adata_sample_col: Column name in adata.obs that contains sample identifiers.

    Returns:
        Updated AnnData object with added sample metadata.
    """
    
    # Store original index
    original_index = adata.obs.index
    
    # Reset index to merge properly
    adata.obs = adata.obs.reset_index(drop=True)

    # Convert categorical columns to string to ensure merge compatibility
    if adata.obs[adata_sample_col].dtype.name == 'category':
        adata.obs[adata_sample_col] = adata.obs[adata_sample_col].astype(str)

    # DEBUG: Check data before merge
    print(f"DEBUG - adata.obs shape before merge: {adata.obs.shape}")
    print(f"DEBUG - sample_metadata shape: {sample_metadata.shape}")
    print(f"DEBUG - adata.obs[{adata_sample_col}] sample values: {adata.obs[adata_sample_col].head(3).tolist()}")
    print(f"DEBUG - sample_metadata[{metadata_sample_col}] sample values: {sample_metadata[metadata_sample_col].head(3).tolist()}")

    # Check for mismatching samples between adata.obs and sample_metadata
    mismatching_samples = set(adata.obs[adata_sample_col]).difference(set(sample_metadata[metadata_sample_col]))
    print(f"DEBUG - Number of mismatching samples: {len(mismatching_samples)}")
    if len(mismatching_samples) > 0:
        print(f"DEBUG - Mismatching samples: {list(mismatching_samples)}")

    # Merge sample metadata with adata.obs based on sample identifiers and region information
    adata.obs = adata.obs.merge(
        sample_metadata,
        left_on=adata_sample_col,
        right_on=metadata_sample_col,
        how="left"
    )
    
    # DEBUG: Check data after merge
    print(f"DEBUG - adata.obs shape after merge: {adata.obs.shape}")

    adata.obs.index = original_index

    return adata

# Main function to process each region and add sample metadata
def main(args: argparse.Namespace):
    
    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    sample_metadata_file = args.sample_metadata
    metadata_region_col = args.metadata_region_col
    metadata_sample_col = args.metadata_sample_col
    adata_sample_col = args.adata_sample_col
    adata_output_prefix = args.adata_output_prefix

    # Load sample metadata    
    print(f"Loading sample metadata from: {sample_metadata_file}")
    sample_metadata = load_sample_metadata(sample_metadata_file, metadata_region_col, metadata_sample_col)
    
    # Iterate over regions and process each one
    for region in regions:
        print(f"Annotating region: {region}")

        # Load region-specific AnnData
        adata_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Loading adata for region {region} from file {adata_fp}")
        adata = sc.read_h5ad(adata_fp)
        print(f"Loaded {adata_fp} with {adata.n_obs} cells and {adata.n_vars} genes")
        
        # Add sample metadata
        adata = add_sample_metadata(adata, sample_metadata, metadata_region_col, metadata_sample_col, adata_sample_col)

        # Save updated AnnData
        output_fp = f"{adata_output_prefix}_{region}.h5ad"
        print(f"Writing filtered adata to {output_fp}")
        adata.write_h5ad(output_fp)
        print(f"Wrote filtered adata to {output_fp}")

        # Clean up memory
        del adata
        gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add sample metadata to region-specific AnnData objects"
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of regions to process",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for input AnnData files (e.g. 'output/region_processed')",
    )
    parser.add_argument(
        "--sample-metadata",
        type=str,
        required=True,
        help="CSV file containing sample metadata with columns for sample identifiers and regions",
    )
    parser.add_argument(
        "--metadata-region-col",
        type=str,
        required=True,
        help="Column name in sample_metadata and adata.obs that contains region information",
    )
    parser.add_argument(
        "--metadata-sample-col",
        type=str,
        required=True,
        help="Column name in sample_metadata that contains sample identifiers",
    )
    parser.add_argument(
        "--adata-sample-col",
        type=str,
        required=True,
        help="Column name in adata.obs that contains sample identifiers (default: 'sample_id')",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Prefix for output AnnData files (e.g. 'output/region_metadata_added')",
    )

    args = parser.parse_args()
    main(args)