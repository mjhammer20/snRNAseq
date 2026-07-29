# Imports

import argparse
import gc
import scanpy as sc
from pathlib import Path
from anndata import AnnData

def gen_qc_plots(
        adata: AnnData,
        qc_plots_dir: str,
        region: str,
        metrics: list
        ) -> None:
    """
    Generate QC plots for the given AnnData object.

    Args:
        adata: The AnnData object to generate QC plots for.
        qc_plots_dir (str): Directory to save the QC plots.
        region (str): The region name for naming the output files.

    Returns:
        None
    """
    
    # Create plots directory if it doesn't exist
    fig_dir = Path(qc_plots_dir)
    if not fig_dir.exists():
        fig_dir.mkdir(parents=True, exist_ok=True)

    # Set CPUs to use for parallel computing
    sc.settings.n_jobs = -1

    sc.settings.verbosity = 1
    sc.settings.figdir = fig_dir
    sc.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=(12, 8)
    )

    # Generate QC plots for the specified metrics
    print(f"Generating QC plots for region: {region}...")
    for metric in metrics:
        print(f"  Plotting metric: {metric}...")
        sc.pl.violin(adata, keys=metric, size=0, save="".join("_" + metric + "_" + region))


def main(args: argparse.Namespace):
    # Extract Arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    qc_plots_dir = args.qc_plots_dir
    metrics = args.metrics.split(", ")

    # Iterate over regions and process each one
    for region in regions:
        print(f"Processing region: {region}...")

        # Load the AnnData object for the specified region
        adata_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Loading AnnData object from file: {adata_fp}...")
        adata = sc.read_h5ad(adata_fp, backed="r")
        print(f"Loaded AnnData object with {adata.n_obs} observations and {adata.n_vars} variables.")

        # Generate QC plots for the current region
        gen_qc_plots(adata, qc_plots_dir, region, metrics)

        # Clean up memory
        del adata
        gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge AnnData objects and plot QC metrics")
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space-separated list of regions to process",
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for input AnnData files",
    )
    parser.add_argument(
        "--qc-plots-dir",
        type=str,
        required=True,
        help="Directory to save QC plots",
    )
    parser.add_argument(
        "--metrics",
        type=str,
        required=True,
        help="Comma-separated list of QC metrics to plot",
    )

    args = parser.parse_args()
    main(args)
