# Imports
import scanpy as sc
import os
from pathlib import Path
import pandas as pd
import numpy as np
from anndata import AnnData
import matplotlib.pyplot as plt
import argparse
import gc
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

# Function to filter AnnData object based on cell type and gene expression
def filter_adata(
        adata: AnnData,
        cell_type_assignment_key: str,
        gene_symbols_key: str,
        min_cells: int
    ) -> AnnData:
    """
    Filter AnnData object based on cell type and gene expression.

    Args:
        adata: AnnData object to filter.
        cell_type_assignment_key: Key in adata.obs for cell type assignments.
        gene_symbols_key: Key in adata.var for gene symbols.
        min_cells: Minimum number of cells required for a cell type to be retained.

    Returns:
        Filtered AnnData object.
    """

    # Create gene symbols column in adata.var if it doesn't exist
    if gene_symbols_key not in adata.var.columns:
        adata.var[gene_symbols_key] = adata.var.index.tolist()

    # Filter cell types with minimum number of observations
    cell_type_counts = adata.obs[cell_type_assignment_key].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= min_cells].index.tolist()

    print(f"Cell types with >= {min_cells} cells: {len(valid_cell_types)}")
    print(f"Excluded cell types: {set(adata.obs[cell_type_assignment_key].unique()) - set(valid_cell_types)}")

    # Subset adata to include only valid cell types
    adata_filtered = adata[adata.obs[cell_type_assignment_key].isin(valid_cell_types)].to_memory()

    print(f"Original adata: {adata.shape}")
    print(f"Filtered adata: {adata_filtered.shape}")

    # Clean up memory
    del adata
    gc.collect()

    return adata_filtered

# Function to determine if genes are expressed in each cell type based on mean and fraction thresholds
def is_expressed(
        adata: AnnData,
        gene_list: list[str],
        cell_type_assignment_key: str,
        condition_key: str,
        mean_threshold: float,
        frac_threshold: float
    ) -> pd.Series:
    """
    Determine if genes are expressed in each cell type based on mean and fraction thresholds.

    Args:
        adata: AnnData object containing the gene expression data.
        gene_list: List of genes to check for expression.
        cell_type_assignment_key: Key in adata.obs for cell type assignments.
        condition_key: Key in adata.obs for condition assignments.
        mean_threshold: Minimum mean expression level to consider a gene expressed.
        frac_threshold: Minimum fraction of cells expressing the gene to consider it expressed. 
    
    Returns:
        expressed: Numpy Series of booleans indicating whether each gene is expressed in each cell type.
    """

    # Convert the sparse matrix to a dense format if necessary and create a DataFrame with gene expression data
    X_dense = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    df = pd.DataFrame(X_dense, index=adata.obs_names, columns=adata.var_names)
    
    # Add cell type assignments and condition IDs to the DataFrame
    df[cell_type_assignment_key] = adata.obs[cell_type_assignment_key].values
    df[condition_key] = adata.obs[condition_key].values

    # Subset to gene list early to avoid operating on the full matrix
    genes = gene_list if gene_list is not None else adata.var_names

    # Check all requested genes exist
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {missing}")

    # Compute mean expression and fraction of cells expressing each gene for each cell type and condition
    mean_expr = df.groupby([cell_type_assignment_key, condition_key])[genes].mean()
    frac_expr = df.groupby([cell_type_assignment_key, condition_key])[genes].apply(lambda g: (g > 0).mean())

    # Determine expressed genes
    expressed = (mean_expr >= mean_threshold) & (frac_expr >= frac_threshold)

    return expressed

# Function to generate a dotplot for the expression of nearest genes across cell types
def generate_dotplot(
        adata: AnnData,
        nearest_genes_in_adata: list[str],
        cell_type_assignment_key: str,
        condition_key: str,
        region: str,
        dotplot_fp: str
    ):
    """
    Generate a dotplot for the expression of nearest genes across cell types.

    Args:
        adata: AnnData object containing the gene expression data.
        nearest_genes_in_adata: List of nearest genes to be plotted.
        cell_type_assignment_key: Key in adata.obs for cell type assignments.
        condition_key: Key in adata.obs for condition assignments.
        region: Name of the brain region for the title of the plot.
        dotplot_fp: File path to save the generated dotplot.

    Returns:
        None. The function saves the dotplot as a PNG file in the specified directory.
    """

    # Initialize dotplot figure
    fig = sc.pl.dotplot(
        adata,
        var_names=nearest_genes_in_adata,
        groupby=[cell_type_assignment_key, condition_key],
        use_raw=False,
        standard_scale='var',
        figsize=(30, 10),
        show=False
    )

    # Get the figure and adjust layout
    fig = plt.gcf()
    fig.suptitle(f'Expression of GWAS Nearest Genes Across Cell Types in {region}', fontsize=20, y=1.02)

    # Save the dotplot
    plt.savefig(dotplot_fp, bbox_inches='tight', dpi=300)

# Function to process a specific brain region
def process_region(
        region: str,
        adata_input_prefix: str,
        nearest_genes: list[str],
        cell_type_assignment_key: str,
        gene_symbols_key: str,
        condition_key: str,
        min_cells: int,
        mean_threshold: float,
        frac_threshold: float,
        results_dir: str,
        expressed_genes_prefix: str,
        dotplots_prefix: str,
        worker_id: int,
        cpus_per_worker: int
    ):
    """
    Process a specific brain region to analyze gene expression.

    Args:
        region: Name of the brain region to process.
        adata_input_prefix: Prefix for the input AnnData file path.
        nearest_genes: List of risk SNP nearest genes from GWAS summary statistics.
        cell_type_assignment_key: Key in adata.obs for cell type assignments.
        gene_symbols_key: Key in adata.var for gene symbols.
        condition_key: Key in adata.obs for condition assignments.
        min_cells: Minimum number of cells for a cell type to be included in the analysis.
        mean_threshold: Minimum mean expression level to consider a gene expressed.
        frac_threshold: Minimum fraction of cells expressing the gene to consider it expressed.
        results_dir: Directory to save the results and plots.
        expressed_genes_prefix: Prefix for the expressed genes output file.
        dotplots_prefix: Prefix for the dotplot output file.
        worker_id: ID of the worker processing this region (for logging purposes).
        cpus_per_worker: Number of CPUs allocated to this worker (for logging purposes).

    Returns:
        None. The function saves the results and plots in the specified directory.
    """

    print(f"Worker {worker_id} - Processing region: {region} (using {cpus_per_worker} CPUs)")

    # Load the AnnData object for the specified region
    adata_input_fp = f"{adata_input_prefix}_{region}.h5ad"
    print(f"Worker {worker_id} - Loading AnnData object for region {region} from {adata_input_fp}...")
    adata = sc.read(adata_input_fp)

    # Load nearest genes for the region
    nearest_genes_in_adata = adata.var_names[adata.var.index.isin(nearest_genes)]

    # Filter the AnnData object to include only cell types with at least 50 observations
    print(f"Worker {worker_id} - Filtering AnnData object for region {region} to include only cell types with at least {min_cells} observations...")
    adata = filter_adata(
        adata=adata,
        cell_type_assignment_key=cell_type_assignment_key,
        gene_symbols_key=gene_symbols_key,
        min_cells=min_cells)

    # Determine expressed genes based on thresholds
    print(f"Worker {worker_id} - Determining expressed genes for region {region} based on mean threshold {mean_threshold} and fraction threshold {frac_threshold}...")
    expressed = is_expressed(
        adata=adata, 
        gene_list=nearest_genes_in_adata,
        cell_type_assignment_key=cell_type_assignment_key,
        condition_key=condition_key,
        mean_threshold=mean_threshold,
        frac_threshold=frac_threshold
    )

    # Save expressed genes and their statistics
    expressed_fp = f"{results_dir}/{expressed_genes_prefix}_{region}.tsv"
    print(f"Worker {worker_id} - Saving expressed genes and their statistics for region {region} to {expressed_fp}...")
    expressed.to_csv(expressed_fp, sep='\t')

    # Generate dotplot for the expression of nearest genes across cell types
    print(f"Worker {worker_id} - Generating dotplot for the expression of nearest genes across cell types in region {region}...")
    dotplot_fp = f"{results_dir}/{dotplots_prefix}_{region}.png"
    generate_dotplot(
        adata=adata,
        nearest_genes_in_adata=nearest_genes_in_adata,
        cell_type_assignment_key=cell_type_assignment_key,
        condition_key=condition_key,
        region=region,
        dotplot_fp=dotplot_fp
    )
    print(f"Worker {worker_id} - Dotplot saved to {dotplot_fp}")

    # Clean up memory
    del adata
    del expressed
    gc.collect()


def main(args: argparse.Namespace):

    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    summary_stats_file = args.gwas_summary_stats
    cell_type_assignment_key = args.cell_type_assignment_key
    gene_symbols_key = args.gene_symbols_key
    condition_key = args.condition_key
    min_cells = args.min_cells
    mean_threshold = args.mean_threshold
    frac_threshold = args.frac_threshold
    results_dir = args.results_dir
    expressed_genes_prefix = args.expressed_genes_prefix
    dotplots_prefix = args.dotplots_prefix

    # Extract risk SNP nearest genes from GWAS summary statistics
    print(f"Loading GWAS summary statistics from {summary_stats_file}")
    gwas_df = pd.read_csv(summary_stats_file, sep='\t')
    nearest_genes = gwas_df['Nearest Gene'].unique()
    nearest_genes = np.append(nearest_genes, 'GBA')
    print(f"Extracted {len(nearest_genes)} unique PD risk SNP nearest genes from GWAS summary statistics")
    
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
                nearest_genes=nearest_genes,
                cell_type_assignment_key=cell_type_assignment_key,
                gene_symbols_key=gene_symbols_key,
                condition_key=condition_key,
                min_cells=min_cells,
                mean_threshold=mean_threshold,
                frac_threshold=frac_threshold,
                results_dir=results_dir,
                expressed_genes_prefix=expressed_genes_prefix,
                dotplots_prefix=dotplots_prefix,
                worker_id=worker_id,
                cpus_per_worker=cpus_per_worker
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
    parser = argparse.ArgumentParser(description="Process brain region gene expression data.")
    parser.add_argument(
        "--regions",
        type=str,
        required=True,
        help="Space separated list of regions to process."
    )
    parser.add_argument(
        "--adata-input-prefix",
        type=str,
        required=True,
        help="Prefix for the input AnnData file path."
    )
    parser.add_argument(
        "--gwas-summary-stats",
        type=str,
        required=True,
        help="Path to the GWAS summary statistics file."
    )
    parser.add_argument(
        "--cell-type-assignment-key",
        type=str,
        required=True,
        help="Key in adata.obs for cell type assignments."
    )
    parser.add_argument(
        "--gene-symbols-key",
        type=str,
        required=True,
        help="Key in adata.var for gene symbols."
    )
    parser.add_argument(
        "--condition-key",
        type=str,
        required=True,
        help="Key in adata.obs for condition assignments."
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        required=True,
        help="Minimum number of cells for a cell type to be included in the analysis."
    )
    parser.add_argument(
        "--mean-threshold",
        type=float,
        required=True,
        help="Minimum mean expression level to consider a gene expressed."
    )
    parser.add_argument(
        "--frac-threshold",
        type=float,
        required=True,
        help="Minimum fraction of cells expressing the gene to consider it expressed."
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        required=True,
        help="Directory to save the results and plots."
    )
    parser.add_argument(
        "--expressed-genes-prefix",
        type=str,
        required=True,
        help="Prefix for the expressed genes output file."
    )
    parser.add_argument(
        "--dotplots-prefix",
        type=str,
        required=True,
        help="Prefix for the dotplot output file."
    )

    args = parser.parse_args()
    main(args)