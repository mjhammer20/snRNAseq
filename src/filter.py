# Imports
import scanpy as sc
import argparse
from anndata import AnnData
import gc

# Filter the merged adata object based on QC metrics and doublet scores
def filter_adata(adata: AnnData, region: str, pct_counts_mt_max: int, doublet_score_max: float, total_counts_limits: list, n_genes_by_counts_limits: list) -> AnnData:
    """
    Filter the merged adata object based on QC metrics and doublet scores.
     - Filters:
        - pct_counts_mt <= 10
        - doublet_score < 0.2
        - total_counts between 100 and 100000
        - n_genes_by_counts between 100 and 10000
     - Uses backed mode to read the data in chunks and apply filters without loading the entire dataset into memory at once
     
    Args:
        adata: AnnData object to filter
        region: The region name for the dataset
        pct_counts_mt_max: Maximum percentage of mitochondrial gene counts allowed per cell
        doublet_score_max: Maximum doublet detection score threshold
        total_counts_limits: Minimum and maximum total UMI (unique molecular identifier) counts per cell
        n_genes_by_counts_limits: Minimum and maximum number of genes detected per cell (genes with at least one count)

    Returns:
        Filtered AnnData object

    """
    
    # Define filters based on QC metrics and doublet scores
    print("Filter 1: pct_counts_mt_max")
    keep_mt = adata.obs['pct_counts_mt'] <= pct_counts_mt_max

    print("Filter 2: doublet_score_max")
    keep_doublet = adata.obs['doublet_score'] < doublet_score_max

    print("Filter 3: total_counts_limits")
    keep_total_counts = (adata.obs['total_counts'] >= total_counts_limits[0]) & (adata.obs['total_counts'] <= total_counts_limits[1])

    print("Filter 4: n_genes_by_counts_limits")
    keep_n_genes_by_counts = (adata.obs['n_genes_by_counts'] >= n_genes_by_counts_limits[0]) & (adata.obs['n_genes_by_counts'] <= n_genes_by_counts_limits[1])

    keep_cells = keep_mt & keep_doublet & keep_total_counts & keep_n_genes_by_counts

    # Create a filterered slice of the input adata object using the boolean mask of cells to keep
    num_keep = keep_cells.sum()
    num_total = len(keep_cells)
    print(f"APPLY FILTER: Keeping {num_keep}/{num_total} cells ({100*num_keep/num_total:.1f}%)")
    adata_filtered = adata[keep_cells]
    print(f"Filtered adata has {adata_filtered.n_obs} cells and {adata_filtered.n_vars} genes")

    return adata_filtered


def main(args: argparse.Namespace):

    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    pct_counts_mt_max = args.pct_counts_mt_max
    doublet_score_max = args.doublet_score_max
    total_counts_limits = args.total_counts_limits
    n_genes_by_counts_limits = args.n_genes_by_counts_limits
    adata_output_prefix = args.adata_output_prefix

    # Disable parallel computing
    sc.settings.n_jobs = 1

    # 1. Filter data
    for region in regions:
        print(f"Filtering region: {region}")

        # Load data
        adata_fp = f"{adata_input_prefix}_{region}.h5ad"
        print(f"Loading adata for region {region} from file {adata_fp} in backed mode")
        adata = sc.read_h5ad(adata_fp, backed='r')
        print(f"Loaded {adata_fp} with {adata.n_obs} cells and {adata.n_vars} genes")

        # Filter adata
        adata_filtered = filter_adata(adata, region, pct_counts_mt_max, doublet_score_max, total_counts_limits, n_genes_by_counts_limits, adata_output_prefix)

        # Write the filtered adata to a new file in compressed format
        output_fp = f"{adata_output_prefix}_{region}.h5ad"
        print(f"Writing filtered adata to {output_fp}")
        adata_filtered.write_h5ad(output_fp)
        print(f"Wrote filtered adata to {output_fp}")

        # Step 4: Close the backed object to release file handles and free up memory
        del adata
        del adata_filtered
        gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter")
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
        help="Input file prefix to read AnnData object from",
    )
    parser.add_argument(
        "--pct-counts-mt-max",
        type=int,
        required=True,
        help="Maximum percentage of mitochondrial gene counts allowed per cell [10]"
    )
    parser.add_argument(
        "--doublet-score-max",
        type=float,
        required=True,
        help="Maximum doublet detection score threshold [0.2]"
    )
    parser.add_argument(
        "--total-counts-limits",
        type=int,
        nargs="+",
        required=True,
        help="Minimum and maximum total UMI (unique molecular identifier) counts per cell [100, 100000]"
    )
    parser.add_argument(
        "--n-genes-by-counts-limits",
        type=int,
        nargs="+",
        required=True,
        help="Minimum and maximum number of genes detected per cell (genes with at least one count) [100, 10000]"
    )
    parser.add_argument(
        "--adata-output-prefix",
        type=str,
        required=True,
        help="Output file prefix to save AnnData object to",
    )

    args = parser.parse_args()
    main(args)
