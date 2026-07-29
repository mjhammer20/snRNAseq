# Imports
import pandas as pd
import argparse
import sys

# Get absolute path to workflows directory in order to import helper functions
workflows_root = 'GP2-Expansion/workflows'
sys.path.insert(0, str(workflows_root))
from src.helpers import normalize_region_names

# Create a list of sample AnnData files based on the metadata and save to a TSV file
def list_sample_adata_files(sample_metadata_file, sample_col, region_col, adata_file_col, adata_dir, sample_adata_suffix):
    # Load sample metadata
    metadata_df = pd.read_csv(sample_metadata_file)
    
    # Create a list of sample AnnData files
    sample_adata_files_df = metadata_df[[sample_col, region_col]].copy()
    sample_adata_files_df[adata_file_col] = sample_adata_files_df[sample_col].apply(lambda x: f"{adata_dir}/{x}.{sample_adata_suffix}.h5ad")

    # Normalize region names in the sample metadata
    sample_adata_files_df[region_col] = normalize_region_names(sample_adata_files_df[region_col])

    return sample_adata_files_df

def main(args: argparse.Namespace):
    # Extract arguments
    sample_metadata_file = args.sample_metadata
    sample_col = args.sample_col
    region_col = args.region_col
    adata_file_col = args.adata_file_col
    adata_dir = args.adata_dir
    sample_adata_suffix = args.sample_adata_suffix
    sample_adata_files = args.output_sample_adata_files
    regions_file = args.output_regions_file
    
    # Generate the list of sample AnnData files
    sample_adata_files_df = list_sample_adata_files(sample_metadata_file, sample_col, region_col, adata_file_col, adata_dir, sample_adata_suffix, sample_adata_files, regions_file)

    # Save sample adata files df to a TSV file
    sample_adata_files_df.to_csv(sample_adata_files, sep='\t', index=False)
    print(f"Sample AnnData files saved to {sample_adata_files}")

    # Save unique regions to a text file
    regions = sample_adata_files_df[region_col].unique()
    with open(regions_file, 'w') as f:
        for region in regions:
            f.write(f"{region}\n")
    print(f"Unique regions saved to {regions_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List sample AnnData files")
    parser.add_argument(
        "--sample-metadata",
        required=True,
        help="Path to the sample metadata CSV file"
    )
    parser.add_argument(
        "--sample-col",
        required=True,
        help="Column name for the sample label in the sample metadata TSV file"
    )
    parser.add_argument(
        "--region-col",
        required=True,
        help="Column name for the region label in the sample metadata TSV file"
    )
    parser.add_argument(
        "--adata-file-col",
        required=True,
        help="Column name for the adata file paths in the output sample adata files TSV"
    )
    parser.add_argument(
        "--adata-dir",
        required=True,
        help="Directory containing the AnnData files"
    )
    parser.add_argument(
        "--sample-adata-suffix",
        required=True,
        help="File suffix for the sample AnnData files"
    )
    parser.add_argument(
        "--output-sample-adata-files",
        required=True,
        help="Path to the output TSV file for sample AnnData files"
    )
    parser.add_argument(
        "--output-regions-file",
        required=True,
        help="Path to the output text file for unique regions"
    )
    args = parser.parse_args()

    main(args)