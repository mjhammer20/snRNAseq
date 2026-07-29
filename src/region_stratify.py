# Imports
import argparse
import scanpy as sc
from anndata import AnnData
import pandas as pd
import os
import gc


# Stratify the merged AnnData object by brain region based on sample metadata
def stratify_by_region(adata: AnnData, region_col: str, adata_input_prefix: str, output_dir: str) -> dict:
    """
    Stratify AnnData object by brain region.
    
    Args:
        adata: Merged AnnData object with sample metadata in .obs
        region_metadata: DataFrame with columns [sample_col, region_col]
        region_col: Column name in adata.obs containing brain region
        output_dir: Directory to save stratified AnnData files
    
    Returns:
        Dictionary of {region: AnnData_subset}

    """
    
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

        region_adata = adata[adata.obs[region_col] == region]
        output_file = f"{output_dir}/{adata_input_prefix}_{region.replace(' ', '_')}.h5ad"
        region_adata.write_h5ad(output_file)
        print(f"Wrote {output_file}")

        # Clear memory
        del region_adata
        gc.collect()
    
    return regions

# Main function to run the script
def main(args: argparse.Namespace):
    # Extract arguments
    adata_input_prefix = args.adata_input_prefix
    datasets = args.datasets.split()
    region_col = args.region_col
    output_dir = args.output_dir

    # Iterate over datasets to stratify by region
    for dataset in datasets:
        print(f"Processing dataset: {dataset}")

        # Load merged AnnData
        adata_fp = f"{adata_input_prefix}_{dataset}.h5ad"
        print(f"Loading {adata_fp}...")
        adata = sc.read_h5ad(adata_fp, backed='r')
        print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        
        # Stratify by region
        regions = stratify_by_region(
            adata_input_prefix=adata,
            region_col=region_col,
            output_dir=output_dir
        )
    
    # Save regions list to file
    regions_output_file = os.path.join(args.output_dir, "regions.txt")
    with open(regions_output_file, "w") as f:
        for region in regions:
            f.write(region.replace(" ", "_") + "\n")
    print(f"Wrote regions list to {regions_output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stratify AnnData by brain region")
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Input file prefix to read merged AnnData object from (without dataset suffix)",
    )
    parser.add_argument(
        "--datasets",
        type=str,
        required=True,
        help="Space-separated list of dataset names to read merged AnnData object for",
    )
    parser.add_argument(
        "--region-col",
        type=str,
        default="Brain Region",
        help="Column name in adata.obs with brain region",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for stratified AnnData files",
    )
    
    args = parser.parse_args()
    main(args)