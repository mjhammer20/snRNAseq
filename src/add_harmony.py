
# Imports
import os
import argparse
import gc
import traceback
import numpy as np
import scanpy as sc
from anndata import AnnData
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to add Harmony integration to AnnData object
def add_harmony(adata: AnnData, batch_key: str, region:str) -> None:
    """
    Adds Harmony integration to AnnData object in place

    Args:
        adata: AnnData object to add Harmony integration to
    
    Returns:
        None, modifies adata in place
    """
    print(f"Running Harmony integration for region: {region}")
    # Check for NaN/Inf in X
    if hasattr(adata.X, 'data'):  # sparse matrix
        has_nan = np.isnan(adata.X.data).any()
        has_inf = np.isinf(adata.X.data).any()
    else:  # dense matrix
        has_nan = np.isnan(adata.X).any()
        has_inf = np.isinf(adata.X).any()
    
    if has_nan or has_inf:
        raise ValueError(f"{region}: AnnData.X contains NaN ({has_nan}) or Inf ({has_inf}) values!")
    
    # Check if latent representation exists and has correct dimensions
    if "X_pca" in adata.obsm:
        print(f"Region {region}: X_pca shape: {adata.obsm['X_pca'].shape}, n_obs: {adata.n_obs}")
        if adata.obsm['X_pca'].shape[0] != adata.n_obs:
            raise ValueError(f"Region {region}: X_pca has {adata.obsm['X_pca'].shape[0]} observations but adata has {adata.n_obs}")
    
    sc.external.pp.harmony_integrate(adata, batch_key, basis="X_pca")

# Function to process a single region
def process_region(
        region: str,
        adata_input_prefix: str,
        adata_output_prefix: str,
        output_metadata_file_prefix: str,
        worker_id: int,
        cpus_per_worker: int,
        batch_key: str
        ) -> None:
    """
    Processes a single region by loading the AnnData object, adding Harmony integration, and saving the results
    Args:q
        region: Region to process
        adata_input_prefix: Prefix for input AnnData file (without region and extension)
        adata_output_prefix: Prefix for output AnnData file (without region and extension)
        output_metadata_file_prefix: Prefix for output metadata file (without region and extension)
    
    Returns:
        None, saves results to disk
    """
    try:
        print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

        # Set CPUs for this worker
        sc.settings.n_jobs = cpus_per_worker

        adata_input = f"{adata_input_prefix}_{region}.h5ad"
        adata_output = f"{adata_output_prefix}_{region}.h5ad"
        output_metadata_file = f"{output_metadata_file_prefix}_{region}.csv"

        # Load Adata object
        print(f"Loading AnnData object: {adata_input}")
        adata = sc.read_h5ad(adata_input)

        # Operations done in place
        add_harmony(adata, batch_key, region)

        # 9. Save the adata
        print(f"Saving AnnData object: {adata_output}")
        adata.write_h5ad(filename=adata_output, compression="gzip")

        # Save metadata
        print(f"Saving metadata file: {output_metadata_file}")
        metatable = adata.obs
        metatable["UMAP_1"] = adata.obsm["X_umap"][:, 0]
        metatable["UMAP_2"] = adata.obsm["X_umap"][:, 1]
        metatable.to_csv(output_metadata_file, index=True)

        print(f"Worker {worker_id} - Region {region} complete!")

        # Clean up memory
        del adata
        gc.collect()
        
    except Exception as e:
        print(f"Worker {worker_id} - ERROR processing region {region}: {e}")
        traceback.print_exc()
        raise e


def main(args: argparse.Namespace):

    # Extract arguments
    batch_key = args.batch_key
    adata_input_prefix = args.adata_input_prefix
    adata_output_prefix = args.adata_output_prefix
    output_metadata_file_prefix = args.output_metadata_file_prefix
    regions = args.regions
    
    # Iterate over regions
    regions = regions.split()
    print(f"Processing {len(regions)} regions: {regions}")

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
                region, adata_input_prefix, adata_output_prefix,
                output_metadata_file_prefix, worker_id, cpus_per_worker, batch_key
            )
            futures[future] = region

        # Wait for tasks as they complete (not in order)
        for future in as_completed(futures):
            region = futures[future]
            try:
                future.result()
                print(f"✓ Region {region} completed")
            except Exception as e:
                print(f"✗ Region {region} failed: {e}")
                traceback.print_exc()

    print("All regions processed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add Harmony integration")
    parser.add_argument(
        "--batch-key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--output-metadata-file-prefix",
        type=str,
        required=True,
        help="Output file to write metadata to",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Regions to process, space separated if multiple",
    )

    args = parser.parse_args()
    main(args)
