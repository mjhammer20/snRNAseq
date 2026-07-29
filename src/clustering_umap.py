# Imports
import argparse
import os
import traceback
import gc
import scanpy as sc
from anndata import AnnData
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to get cluster and UMAP
def get_cluster_umap(
        adata: AnnData, 
        latent_key: str, 
        n_neighbors : int, 
        leiden_res: list[float]
    ) -> AnnData:
    """
    Get cluster and UMAP embeddings for the given AnnData object.

    Args:
        adata: AnnData object to get cluster and UMAP for
        latent_key: Key in adata.obsm to use for neighbors calculation and Leiden clustering
        n_neighbors: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation [15]

    Returns:
        AnnData object with Leiden clusters and UMAP embeddings added
    """

    # Calculate neighbor graph on scVI latent
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=latent_key)
    # Do Leiden
    for resolution in leiden_res:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=f"leiden_res_{resolution:.2f}",
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    sc.tl.umap(adata)
    return adata

def process_region(
        region: str, 
        adata_input_prefix: str, 
        adata_output_prefix: str, 
        latent_key: str, 
        n_neighbors: int, 
        leiden_res: list[float], 
        worker_id: int,
        cpus_per_worker: int
    ):
    """
    Process a single region by loading the corresponding AnnData object, getting cluster and UMAP, and saving the processed AnnData object.

    Args:
        region: Region to process
        adata_input_prefix: Prefix for input AnnData file
        adata_output_prefix: Prefix for output AnnData file
        latent_key: Key in adata.obsm to use for neighbors calculation and Leiden clustering
        n_neighbors: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation [15]
        leiden_res: List of Leiden resolutions which are the parameter values controlling the coarseness of the clustering

    Returns:
        None
    """
    try:
        print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

        # Set CPUs for this worker
        sc.settings.n_jobs = cpus_per_worker

        # Load data
        adata_input = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Worker {worker_id} - Loading data from {adata_input}")
        adata = sc.read_h5ad(adata_input)
        
        # Get cluster and UMAP
        print(f"Worker {worker_id} - Generating clusters and UMAP for region {region}")
        adata = get_cluster_umap(adata=adata, latent_key=latent_key, n_neighbors=n_neighbors, leiden_res=leiden_res)
        
        # Save adata
        adata_output = f"{adata_output_prefix}_{region}.h5ad"
        print(f"Worker {worker_id} - Saving processed data to {adata_output}")
        adata.write_h5ad(filename=adata_output, compression="gzip")

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
    latent_key = args.latent_key
    n_neighbors = args.n_neighbors
    leiden_res = args.leiden_res
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    adata_output_prefix = args.adata_output_prefix

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
                region=region, adata_input_prefix=adata_input_prefix, adata_output_prefix=adata_output_prefix,
                latent_key=latent_key, n_neighbors=n_neighbors, leiden_res=leiden_res, worker_id=worker_id, cpus_per_worker=cpus_per_worker
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
    parser = argparse.ArgumentParser(description="Annotate clusters")
    parser.add_argument(
        "--latent-key",
        type=str,
        required=True,
        help="Latent key to save the scVI latent to",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        required=True,
        help="The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation",
    )
    parser.add_argument(
        "--leiden-res",
        type=float,
        nargs="+",
        required=True,
        help="Leiden resolutions which are the parameter values controlling the coarseness of the clustering",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for input AnnData files",
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Prefix for output AnnData files",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Regions to process, separated by space",
    )

    args = parser.parse_args()
    main(args)
