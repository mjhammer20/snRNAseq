# Imports
import argparse
import scanpy as sc
from anndata import concat as ad_concat
from anndata import AnnData
import gc
import pandas as pd
import os

def merge_dataset_sample_adatas(
        sample_adata_files: str,
        dataset_label: str,
        sample_label: str,
        adata_file_label: str,
        dataset_name: str
        ) -> AnnData:
    '''
    Merge all sample AnnData objects for a specific dataset.
    
    Args:
        sample_adata_files (str): Path to the file listing sample adata files (tab-separated).
        dataset_label (str): The column name in the sample_adata_files that contains the dataset labels.
        sample_label (str): The column name in the sample_adata_files that contains the sample labels.
        adata_file_label (str): The column name in the sample_adata_files that contains the file paths for the adata objects.
        dataset_name (str): The specific dataset name to filter samples by.
    
    Returns:
        AnnData: Merged AnnData object for the specified dataset.
    '''

    # Load the sample adata files information
    sample_adata_files_df = pd.read_csv(sample_adata_files, sep='\t')

    # Filter the sample files for the specified dataset
    dataset_sample_files_df = sample_adata_files_df[sample_adata_files_df[dataset_label] == dataset_name]

    # Merge the sample adata files for this dataset
    print(f"Merging dataset '{dataset_name}' with {len(dataset_sample_files_df)} samples")
    missing_files = []

    # Iterate over the sample files for this dataset and merge them one by one to avoid memory issues
    for idx, row in dataset_sample_files_df.iterrows():
        sample_name = row[sample_label]
        file_path = row[adata_file_label]
        print(f"  Loading sample {sample_name} from file {file_path}")
        
        # Check if the file exists before attempting to load it
        if os.path.exists(file_path):

            # For the first sample, initialize the merged adata. For subsequent samples, concatenate with the existing merged adata.
            if idx == dataset_sample_files_df.index[0]:  # First sample initializes the merged adata
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
        print(f"Missing files for dataset '{dataset_name}': {missing_files}")

    return merged_adata


def main(args: argparse.Namespace):
    # Extract arguments
    sample_adata_files = args.sample_adata_files
    adata_output_prefix = args.adata_output_prefix
    output_metadata_prefix = args.output_metadata_prefix
    dataset_label = args.dataset_label
    sample_label = args.sample_label
    adata_file_label = args.adata_file_label
    datasets = args.datasets

    # Parse space-separated datasets string
    datasets = datasets.split()

    # Iterate over datasets and merge sample adata files
    for dataset in datasets:

        # Define output file paths for merged adata and metadata
        adata_output = f'{adata_output_prefix}_{dataset}.h5ad'
        metadata_output = f'{output_metadata_prefix}_{dataset}.csv'

        # Merge adatas
        merged_adata = merge_dataset_sample_adatas(sample_adata_files, dataset_label, sample_label, adata_file_label, dataset)

        # Save merged adata
        print(f"Saving merged AnnData object: {adata_output}")
        merged_adata.write_h5ad(filename=adata_output, compression="gzip")  

        # Save metadata
        print(f"Saving metadata file: {metadata_output}")
        metadata = merged_adata.obs
        metadata.to_csv(metadata_output, index=True)

        # Clean up memory
        del merged_adata
        gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge sample adata files in chunks and save merged adata and metadata.")
    parser.add_argument(
        "--sample-adata-files",
        type=str,
        required=True,
        help="Path to the file listing sample adata files (tab-separated)."
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        default="merged_adata_chunk",
        help="Prefix for output merged AnnData files."
    )
    parser.add_argument(
        "--output-metadata-prefix",
        type=str,
        required=True,
        help="Prefix for output metadata CSV files."
    )
    parser.add_argument(
        "--dataset-label",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the dataset labels."
    )
    parser.add_argument(
        "--sample-label",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the sample labels."
    )
    parser.add_argument(
        "--adata-file-label",
        type=str,
        required=True,
        help="Column name in the sample adata files that contains the file paths for the adata objects."
    )
    parser.add_argument(
        "--datasets",
        type=str,
        required=True,
        help="Space-separated list of dataset names to process (e.g., 'dataset1 dataset2 dataset3')."
    )
    
    args = parser.parse_args()
    main(args)

