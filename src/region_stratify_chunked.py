# Imports
import argparse
import scanpy as sc
from anndata import AnnData
import pandas as pd
import os
import gc

# Load metadata mapping samples to brain regions
def load_region_metadata(file_path: str) -> pd.DataFrame:
    """
    Load region metadata from a comma-separated values (CSV) file.

    Args:
        file_path: Path to the CSV file containing region metadata.

    Returns:
        DataFrame containing region metadata.
    """
    region_metadata = pd.read_csv(file_path)
    return region_metadata

# Normalize region names to eliminate inconsistencies in metadata
def normalize_region_names(region_series: pd.Series) -> pd.Series:
    """
    Normalize region names by stripping whitespace and converting to title case.
    
    Args:
        region_series: pandas Series with region names
    
    Returns:
        Normalized region names
    """
    return region_series.str.strip().str.title()

# Stratify the merged AnnData object by brain region based on sample metadata
def stratify_by_region(adata: AnnData, region_metadata: pd.DataFrame, region_sample_col: str, adata_sample_col: str, region_col: str, output_dir: str, output_suffix: str ) -> None:
    """
    Stratify AnnData object by brain region.
    
    Args:
        adata: Merged AnnData object with sample metadata in .obs
        region_metadata: DataFrame with columns [sample_col, region_col]
        sample_col: Column name in adata.obs containing sample identifiers
        region_col: Column name in region_metadata containing brain region
        output_dir: Directory to save stratified AnnData files
        output_suffix: Suffix for output file names
    
    Returns:
        None

    """

    # Normalize region names in the metadata
    region_metadata[region_col] = normalize_region_names(region_metadata[region_col])
    
    # Add region information to adata.obs
    print(f"Adding region information to adata.obs...")
    adata.obs = adata.obs.reset_index(drop=True)
    adata.obs = adata.obs.merge(
        region_metadata[[region_sample_col, region_col]],
        left_on=adata_sample_col,
        right_on=region_sample_col,
        how="left"
    )
    
    # Write metadata file for reference
    metadata_output_file = os.path.join(output_dir, f"merged_cleaned_filtered_{output_suffix}_metadata.csv")
    adata.obs.to_csv(metadata_output_file)
    print(f"  Saved: {metadata_output_file}")

    # Check for samples without region assignment
    missing_regions = adata.obs[adata.obs[region_col].isna()]
    if len(missing_regions) > 0:
        print(f"Warning: {len(missing_regions)} cells from samples without region assignment")
    
    # Stratify by region
    print(f"Stratifying by brain region...")
    regions = adata.obs[region_col].unique()
    regions = [r for r in regions if pd.notna(r)]  # remove NaN
    
    for region in sorted(regions):

        # Filter to region
        region_adata = adata[adata.obs[region_col] == region]
        print(f"  {region}: {region_adata.n_obs} cells, {region_adata.n_vars} genes")

        # Save region-specific AnnData
        region_dir = os.path.join(output_dir, region.replace(" ", "_"))
        os.makedirs(region_dir, exist_ok=True)
        output_file = os.path.join(region_dir, f"merged_cleaned_filtered_{output_suffix}.h5ad")
        region_adata.write_h5ad(output_file, compression="gzip")
        print(f"  Saved: {output_file}")

        # Clean up memory
        del region_adata
        gc.collect()
    

# Main function to run the script
def main(args: argparse.Namespace):

    # Load region metadata
    print(f"Loading region metadata from {args.region_metadata}...")
    region_metadata = load_region_metadata(args.region_metadata)
    print(f"Region metadata columns: {list(region_metadata.columns)}")
    
    # Load merged AnnData objects (assuming they are split into chunks)
    suffixes = ['12', '22', '32']

    for suffix in suffixes:
        adata_path = f'{args.output_dir}/{args.adata_prefix}_{suffix}.h5ad'
        # Load merged AnnData
        print(f"Loading {adata_path}...")
        adata = sc.read_h5ad(adata_path)
        print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    
        # Stratify by region
        stratify_by_region(
            adata = adata,
            region_metadata = region_metadata,
            region_sample_col = args.region_sample_col,
            adata_sample_col = args.adata_sample_col,
            region_col = args.region_col,
            output_dir = args.output_dir,
            output_suffix = suffix
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stratify AnnData by brain region")
    parser.add_argument(
        "--adata-prefix",
        type=str,
        required=True,
        help="Input merged AnnData file",
    )
    parser.add_argument(
        "--region-metadata",
        type=str,
        required=True,
        help="CSV file with sample-to-region mapping",
    )
    parser.add_argument(
        "--region-sample-col",
        type=str,
        default="Asap Sample ID",
        help="Column name in adata.obs with sample identifiers",
    )
    parser.add_argument(
        "--adata-sample-col",
        type=str,
        default="sample",
        help="Column name in region metadata with sample identifiers",
    )
    parser.add_argument(
        "--region-col",
        type=str,
        default="Brain Region",
        help="Column name in region metadata with brain region",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for stratified AnnData files",
    )
    args = parser.parse_args()
    main(args)