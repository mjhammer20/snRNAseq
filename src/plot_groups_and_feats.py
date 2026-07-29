# Imports
import argparse
import os 
import gc 
import traceback
from pathlib import Path
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to make feature and group UMAP plots
def make_f_and_g_plots(adata: sc.AnnData, groups:str, features:str, region:str):
    """
    Plot UMAPs colored by groups and features.

    Args:
        adata: AnnData object containing the data to plot
        groups: Comma-separated list of group keys in adata.obs to plot
        features: Comma-separated list of feature keys in adata.obs to plot
        region: Region being processed

    Returns:
        None; saves UMAP plots to disk
    """

    features_list = features.split(", ")
    plot_features = [x for x in features_list if x in adata.obs.columns]
    file_name = f"_features_{region}.png"
    sc.pl.embedding(
        adata,
        basis="umap",
        color=plot_features,
        frameon=False,
        show=False,
        ncols=1,
        save=file_name,
    )

    groups_list = groups.split(", ")
    plot_groups = [x for x in groups_list if x in adata.obs.columns]
    file_name = f"_groups_{region}.png"
    sc.pl.embedding(
        adata,
        basis="umap",
        color=plot_groups,
        frameon=False,
        show=False,
        ncols=1,
        save=file_name,
    )

def process_region(adata_input_prefix: str, region: str, groups:str, features:str, worker_id: int, cpus_per_worker: int, umap_dir: str):
    """
    Process a single region by loading the corresponding AnnData object, plotting UMAPs colored by groups and features, and saving the results

    Args:
        adata_input_prefix: Prefix for input AnnData file
        region: Region to process
        groups: Comma-separated list of group keys in adata.obs to plot
        features: Comma-separated list of feature keys in adata.obs to plot
        worker_id: ID of the worker process (for logging purposes)
        cpus_per_worker: Number of CPUs to use for this worker (for logging purposes)
        umap_dir: Output directory to save UMAP plots to

    Returns:
        None; saves UMAP plots to disk
    """
    print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

    # Load data
    adata_input = f"{adata_input_prefix}_{region}.h5ad"
    print(f"Loading AnnData object: {adata_input}")
    adata = sc.read_h5ad(adata_input)

    # Set CPUs for this worker
    sc.settings.n_jobs = cpus_per_worker

    # Create report directory if it doesn't exist
    if not Path(umap_dir).exists():
        Path(umap_dir).mkdir(parents=True, exist_ok=True)

    # Set working directory and load packages
    sc.settings.verbosity = 1
    sc.settings.figdir = umap_dir
    sc.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
    )  # type: ignore

    # Make plots
    print(f"Making UMAP plots for region: {region}")
    make_f_and_g_plots(adata=adata, groups=groups, features=features, region=region)

    print(f"Worker {worker_id} - Region {region} complete!")

    # Clean up memory
    del adata
    gc.collect()


def main(args: argparse.Namespace):
    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    groups = args.groups
    features = args.features
    umap_dir = args.umap_dir

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
                adata_input_prefix=adata_input_prefix, region=region, groups=groups, features=features,
                worker_id=worker_id, cpus_per_worker=cpus_per_worker, umap_dir=umap_dir
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
    parser = argparse.ArgumentParser(description="Plot groups and features")
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for input AnnData file",
    )
    parser.add_argument(
        "--groups",
        type=str,
        required=True,
        help="Group to plot umaps for"
    )
    parser.add_argument(
        "--features",
        type=str,
        required=True,
        help="Feature to plot umaps for"
    )
    parser.add_argument(
        "--umap-dir",
        type=str,
        required=True,
        help="Output directory to save UMAP plots to",
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Regions to process, space separated if multiple",
    )

    args = parser.parse_args()
    main(args)
