# Imports
import os
import traceback
import argparse
import gc
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation
from pathlib import Path
from anndata import AnnData
from pandas import DataFrame
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to compute `scib` metrics on final artifacts
def get_artifact_metrics(
    adata: AnnData, 
    batch_key: str, 
    cell_type_assignment_key: str, 
    scib_report_dir: Path | str
) -> DataFrame:
    """
    Compute `scib` metrics on final artifacts

    Args:
        adata: AnnData object containing the final artifacts to evaluate with `scib`
        batch_key: Key in adata.obs for batch information
        cell_type_assignment_key: Key in adata.obs for cell type assignments
        scib_report_dir: Directory to save `scib` report
    Returns:
        DataFrame containing `scib` metrics results
    """

    # Check which embeddings exist
    print(f"Available embeddings in obsm: {list(adata.obsm.keys())}")
    
    # Only use embeddings that exist
    embedding_obsm_keys = []

    if "X_umap" in adata.obsm:
        embedding_obsm_keys.append("X_umap")

    if not embedding_obsm_keys:
        raise ValueError("No embeddings found in adata.obsm!")
    
    print(f"Using embeddings: {embedding_obsm_keys}")

    biocons = BioConservation(isolated_labels=False)
    
    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=cell_type_assignment_key,
        embedding_obsm_keys=embedding_obsm_keys,
        bio_conservation_metrics=biocons,
        batch_correction_metrics=None,  # Skip batch correction metrics to avoid concatenation errors
        n_jobs=-1,
    )
    bm.benchmark()

    bm.plot_results_table(min_max_scale=False, save_dir=scib_report_dir, )
    df = bm.get_results(min_max_scale=False)

    return df


def process_region(
        region: str,
        adata_input_prefix: str,
        batch_key: str,
        cell_type_assignment_key: str,
        scib_report_dir: str, 
        worker_id: int,
        cpus_per_worker: int
) -> None:
    """
    Process a single region by loading the corresponding AnnData object, computing `scib` metrics, and saving the results

    Args:
        region: Region to process
        adata_input_prefix: Prefix for input AnnData file
        batch_key: Key in AnnData object for batch information
        cell_type_assignment_key: Key in AnnData object for cell type assignments
        scib_report_dir: Prefix for output directory to save `scib` report
    Returns:
        None, saves results to disk
    """

    print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

    # Set CPUs for this worker
    sc.settings.n_jobs = cpus_per_worker
    
    # Load data
    adata_input = f"{adata_input_prefix}_{region}.h5ad"
    print(f"Loading AnnData object: {adata_input}")
    adata = sc.read_h5ad(adata_input)  # type: ignore

    # Create report directory if it doesn't exist
    region_report_dir = f"{scib_report_dir}/{region}"
    if not Path(region_report_dir).exists():
        Path(region_report_dir).mkdir(parents=True, exist_ok=True)

    # Compute `scib` metrics and save report
    print(f"Computing `scib` metrics for region: {region}")
    report_df = get_artifact_metrics(adata=adata, batch_key=batch_key, cell_type_assignment_key=cell_type_assignment_key, scib_report_dir=region_report_dir)
    report_df_fp = f"{region_report_dir}/scib_report.csv"
    report_df.to_csv(report_df_fp, index=True)
    print(f"Worker {worker_id} - Region {region} complete! Report saved to: {report_df_fp}")

    # Clean up memory
    del adata
    gc.collect()


def main(args: argparse.Namespace):
    # Extract arguments
    regions = args.regions.split()
    cell_type_assignment_key = args.cell_type_assignment_key    
    batch_key = args.batch_key
    adata_input_prefix = args.adata_input_prefix
    output_report_dir = args.output_report_dir

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
                region=region, adata_input_prefix=adata_input_prefix, batch_key=batch_key, cell_type_assignment_key=cell_type_assignment_key,
                scib_report_dir=output_report_dir, worker_id=worker_id, cpus_per_worker=cpus_per_worker
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
    parser = argparse.ArgumentParser(
        description="Compute `scib` metrics on final artifacts"
    )
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of regions to process."
    )
    parser.add_argument(
        "--cell-type-assignment-key",
        type=str,
        required=True,
        help="scANVI cell type predictions column name in AnnData object that references our 'cell_type' labels",
    )
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
        help="Prefix for input AnnData files (without region suffix)",
    )
    parser.add_argument(
        "--output-report-dir",
        type=str,
        required=True,
        help="Output folder to save `scib` report",
    )

    args = parser.parse_args()
    main(args)
