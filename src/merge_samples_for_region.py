# Imports
import argparse
import scanpy as sc
from anndata import concat as ad_concat
from anndata import AnnData
import gc
import pandas as pd
import os

# Function to merge sample adata files for a single region
def merge_region_sample_adatas(
        sample_adata_files: str,
        region_col: str,
        sample_col: str,
        adata_file_col: str,
        region: str
        ) -> AnnData:
    '''
    Merge all sample AnnData objects for a specific region.

    Args:
        sample_adata_files (str): Path to the file listing sample adata files (tab-separated).
        region_col (str): The column name in the sample_adata_files that contains the dataset labels.
        sample_col (str): The column name in the sample_adata_files that contains the sample labels.
        adata_file_col (str): The column name in the sample_adata_files that contains the file paths for the adata objects.
        region (str): The specific region to filter samples by.
    
    Returns:
        AnnData: Merged AnnData object for the specified region.
    '''

    # Load the sample adata files information
    sample_adata_files_df = pd.read_csv(sample_adata_files, sep='\t')

    # Filter the sample files for the specified region
    region_sample_files_df = sample_adata_files_df[sample_adata_files_df[region_col] == region]

    # Merge the sample adata files for this region
    print(f"Merging region '{region}' with {len(region_sample_files_df)} samples")
    missing_files = []

    # Iterate over the sample files for this region and merge them one by one to avoid memory issues
    for idx, row in region_sample_files_df.iterrows():
        sample_name = row[sample_col]
        file_path = row[adata_file_col]
        print(f"  Loading sample {sample_name} from file {file_path}")
        
        # Check if the file exists before attempting to load it
        if os.path.exists(file_path):

            # For the first sample, initialize the merged adata. For subsequent samples, concatenate with the existing merged adata.
            if idx == region_sample_files_df.index[0]:  # First sample initializes the merged adata
                merged_adata = sc.read_h5ad(file_path)
            else:
                adata = sc.read_h5ad(file_path)
                merged_adata = ad_concat([merged_adata, adata], merge="same", uns_merge="same", index_unique="_")
            
                # Clean up memory
                del adata
                gc.collect()
        else:
            missing_files.append(file_path)

    # Log any missing files
    if missing_files:
        print(f"Missing files for dataset '{region}': {missing_files}")

    return merged_adata

# Main function to process each region and save merged adata and metadata
def main(args: argparse.Namespace):
    # Extract arguments
    regions = args.regions.split()
    sample_adata_files = args.sample_adata_files
    region_col = args.region_col
    sample_col = args.sample_col
    adata_file_col = args.adata_file_col
    adata_output_prefix = args.adata_output_prefix
    output_metadata_prefix = args.output_metadata_prefix
    
    # Iterate over regions and process each one
    for region in regions:
        print(f"Processing region: {region}")
        
        # Define output file paths for merged adata and metadata
        adata_output_fp = f'{adata_output_prefix}_{region}.h5ad'
        metadata_output_fp = f'{output_metadata_prefix}_{region}.csv'

        # Merge adatas for the single region
        merged_adata = merge_region_sample_adatas(sample_adata_files, region_col, sample_col, adata_file_col, region)

        # Save merged adata
        print(f"Writing merged AnnData object: {adata_output_fp}")
        merged_adata.write_h5ad(filename=adata_output_fp, compression="gzip")
        print(f"Saved merged AnnData object: {adata_output_fp}")

        # Save metadata (obs) if present; otherwise write an empty csv with header
        print(f"Writing metadata file: {metadata_output_fp}")
        merged_adata.obs.to_csv(metadata_output_fp, index=True)
        print(f"Saved metadata file: {metadata_output_fp}")

        # Clean up memory
        del merged_adata
        gc.collect()

        print(f"Finished processing region: {region}")
    
    print("All regions processed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge sample adata files for a single region and save merged adata and metadata.")
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of region names to process."
    )
    parser.add_argument(
        "--sample-adata-files",
        type=str,
        required=True,
        help="Path to the file listing sample adata files (tab-separated)."
    )
    parser.add_argument(
        "--region-col",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the region labels."
    )
    parser.add_argument(
        "--sample-col",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the sample labels."
    )
    parser.add_argument(
        "--adata-file-col",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the file paths for the adata objects."
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        default="merged_adata_chunk",
        help="Prefix for output merged AnnData files. The dataset name will be appended."
    )
    parser.add_argument(
        "--output-metadata-prefix",
        type=str,
        required=True,
        help="Prefix for output metadata CSV files. The dataset name will be appended."
    )
    
    args = parser.parse_args()
    main(args)
