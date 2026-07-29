# Imports
import argparse
import scanpy as sc
from anndata import AnnData
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc
import os

# Function to normalize and log-transform the gene expression data in the AnnData object
def norm_transform(adata: AnnData, target_sum: float) -> AnnData:
    """
    Normalize and log-transform the gene expression data in the AnnData object.

    Args:
        adata: AnnData object to process
        target_sum: Target sum for normalization (e.g., 1e6 for CPM)

    Returns:
        AnnData: Processed AnnData object with normalized and log-transformed data

    """
    # Retain raw counts in a separate layer
    adata.layers["counts"] = adata.X.copy()

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    return adata

# Function to add cell type assignment labels to the AnnData object for each region
def assign_cell_type_assignment_labels(
        adata: AnnData, 
        mmc_cell_type_label_key: str,
        scanvi_prediction_key: str,
        cell_type_assignment_key: str
        ) -> AnnData:
    """
    Add cell type assignment labels to the AnnData object for each region

    Args:
        adata: AnnData object to add cell type labels to
        mmc_cell_type_label_key: Key in adata.obs to use for cell type labels
        scanvi_prediction_key: Key in adata.obs to use for scANVI predictions
        cell_type_assignment_key: Key in adata.obs to save the final cell type assignment labels

    Returns:
        AnnData object with cell type assignment labels added
    """

    # Assign cell type labels based on scANVI predictions and MMC results
    adata.obs[cell_type_assignment_key] = adata.obs.apply(
        lambda row: row[mmc_cell_type_label_key] if row[mmc_cell_type_label_key] != "Unknown" else row[scanvi_prediction_key], axis=1)

    return adata

def process_region(
        region: str, 
        adata_input_prefix: str,
        target_sum: float,
        mmc_cell_type_label_key: str,
        scanvi_prediction_key: str,
        cell_type_assignment_key: str,
        adata_output_prefix: str,
        metadata_output_prefix: str,
        worker_id: int,
        cpus_per_worker: int
    ):
    """
    Process a single region by loading the corresponding AnnData object, assigning cell type labels, and saving the processed AnnData object.

    Args:
        region: Region to process
        adata_input_prefix: Prefix for input AnnData file
        target_sum: Target sum for normalization (e.g., 1e6 for CPM)
        mmc_cell_type_label_key: Key in adata.obs to use for cell type labels
        scanvi_prediction_key: Key in adata.obs to use for scANVI predictions
        cell_type_assignment_key: Key in adata.obs to save the final cell type assignment labels
        adata_output_prefix: Prefix for output AnnData file
        metadata_output_prefix: Prefix for output metadata file
        worker_id: ID of the worker processing this region
        cpus_per_worker: Number of CPUs allocated to this worker
    """

    print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

    # Load the AnnData object
    adata_input_fp = f"{adata_input_prefix}_{region}.h5ad"
    print(f"Worker {worker_id} - Loading AnnData object for region {region} from {adata_input_fp}...")
    adata = sc.read_h5ad(adata_input_fp)

    # Normalize and log-transform the data
    print(f"Worker {worker_id} - Normalizing and log-transforming data for region {region}...")
    adata = norm_transform(adata, target_sum=target_sum)

    # Assign cell type labels
    print(f"Worker {worker_id} - Assigning cell type labels for region {region}...")
    adata = assign_cell_type_assignment_labels(
        adata=adata,
        mmc_cell_type_label_key=mmc_cell_type_label_key,
        scanvi_prediction_key=scanvi_prediction_key,
        cell_type_assignment_key=cell_type_assignment_key
    )

    # Save the processed AnnData object
    adata_output_fp = f"{adata_output_prefix}_{region}.h5ad"
    print(f"Worker {worker_id} - Saving AnnData object with assigned cell types for region {region} to {adata_output_fp}...")
    adata.write_h5ad(adata_output_fp, compression="gzip")

    # Save adata.obs to a TSV file
    metadata_output_fp = f"{metadata_output_prefix}_{region}.tsv"
    print(f"Worker {worker_id} - Saving metadata for region {region} to {metadata_output_fp}...")
    adata.obs.to_csv(metadata_output_fp, sep="\t", index=False)

    # Clean up memory
    del adata
    gc.collect()

def main(args: argparse.Namespace):
    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    target_sum = args.target_sum
    mmc_cell_type_label_key = args.mmc_cell_type_label_key
    scanvi_prediction_key = args.scanvi_prediction_key
    cell_type_assignment_key = args.cell_type_assignment_key
    adata_output_prefix = args.adata_output_prefix
    metadata_output_prefix = args.metadata_output_prefix

    # Number of parallel workers
    n_workers = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))

    # Allocate CPUs
    total_cpus = os.cpu_count()
    cpus_per_worker = max(1, total_cpus // n_workers)
    
    print(f"Total CPUs: {total_cpus}, CPUs per worker: {cpus_per_worker}")

    # Parallelize across regions
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {}
        for idx, region in enumerate(regions):
            worker_id = idx % n_workers
            future = executor.submit(
                process_region,
                region=region,
                adata_input_prefix=adata_input_prefix,
                target_sum=target_sum,
                mmc_cell_type_label_key=mmc_cell_type_label_key,
                scanvi_prediction_key=scanvi_prediction_key,
                cell_type_assignment_key=cell_type_assignment_key,
                adata_output_prefix=adata_output_prefix,
                metadata_output_prefix=metadata_output_prefix,
                worker_id=worker_id,
                cpus_per_worker=cpus_per_worker
            )
            futures[future] = region

        for future in as_completed(futures):
            region = futures[future]
            try:
                future.result()
                print(f"Successfully processed region {region}")
            except Exception as e:
                print(f"Error processing region {region}: {e}")

        print("All regions processed!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign cell type labels to AnnData objects for each region"
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of regions to process"
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for input AnnData files (e.g., 'adata_region')"
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        required=True,
        help="Target sum for normalization (e.g., 1e6 for CPM)"
    )
    parser.add_argument(
        "--mmc-cell-type-label-key",
        type=str,
        required=True,
        help="Key in adata.obs to use for MMC cell type labels"
    )
    parser.add_argument(
        "--scanvi-prediction-key",
        type=str,
        required=True,
        help="Key in adata.obs to use for scANVI predictions"
    )
    parser.add_argument(
        "--cell-type-assignment-key",
        type=str,
        required=True,
        help="Key in adata.obs to save the final cell type assignment labels"
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Prefix for output AnnData files (e.g., 'adata_region_with_labels')"
    )
    parser.add_argument(
        "--metadata-output-prefix",
        type=str,
        required=True,
        help="Prefix for output metadata files"
    )
    
    args = parser.parse_args()
    main(args)