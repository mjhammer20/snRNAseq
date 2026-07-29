# Imports
import gc
import os
import traceback
import pandas as pd
import argparse
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from concurrent.futures import ProcessPoolExecutor, as_completed


def initialize_dds(adata, inference):
    """
    Initialize DESeq2 dataset from AnnData object.
    Args:
        adata (AnnData): Annotated data object.
        inference (DefaultInference): DESeq2 inference object.
    Returns:
        dds (DESeqDataSet): DESeq2 dataset object.
    """
    # Extract counts and metadata
    counts_array = adata.layers["counts"]
    if hasattr(counts_array, 'toarray'):
        counts_array = counts_array.toarray()
    metadata = adata.obs.copy()
    
    # Convert counts to DataFrame (cells x genes)
    counts_df = pd.DataFrame(
        counts_array,
        index=adata.obs_names,
        columns=adata.var_names
    )
    
    # Pseudobulk by summing counts across cells for each sample
    pb_counts = counts_df.groupby(by=metadata["sample"], axis=0).sum()  # samples x genes

    # Get unique metadata per sample using agg with 'first' for all columns
    pb_metadata = metadata.groupby("sample", observed=True).agg({col: 'first' for col in metadata.columns})

    # Convert CategoricalIndex to regular Index to remove empty categories
    pb_metadata.index = pd.Index(pb_metadata.index.astype(str))

    # Verify alignment before proceeding
    if pb_counts.shape[0] != pb_metadata.shape[0]:
        print(f"ERROR: Shape mismatch! pb_counts: {pb_counts.shape[0]}, pb_metadata: {pb_metadata.shape[0]}")
        raise ValueError("pb_counts and pb_metadata row counts don't match")

    # Initialize DESeq2 dataset
    dds = DeseqDataSet(
        counts=pb_counts,
        metadata=pb_metadata,
        design_factors=["`Condition ID`"],
        ref_level=["`Condition ID`", "Control"],
        inference=inference, 
        refit_cooks=True
    )

    return dds

def run_deseq2(dds, inference):
    """
    Run DESeq2 analysis on dataset.
    Args:
        dds (DESeqDataSet): DESeq2 dataset object.
        inference (DefaultInference): DESeq2 inference object.
    Returns:
        results (DataFrame): DESeq2 results table.
    """
    # Fit dispersions and LFC
    dds.deseq2()

    # Get results
    ds = DeseqStats(dds, contrast=["Condition ID", "PD", "Control"], inference=inference)
    ds.summary()
    results = ds.results_df

    return results

def identify_significant_genes(results):
    """
    Identify significant genes based on adjusted p-value and log2 fold change thresholds.
    Args:
        results (DataFrame): DESeq2 results table.
    Returns:
        significant_genes (DataFrame): Subset of results with significant genes.
    """
    # Filter for significance
    significant_genes = results[    
        (results["padj"] < 0.05) & 
        (results["log2FoldChange"].abs() > 1)
    ]
    return significant_genes

def process_region(region, adata_input_prefix, results_prefix, cpus_per_worker):
    """
    Process a single region by loading the corresponding AnnData object, running DESeq2 analysis, and saving the results.
    Args:
        region (str): Region to process.
        adata_input_prefix (str): Prefix for input AnnData file.
        results_prefix (str): Prefix for output results file.
        cpus_per_worker (int): Number of CPUs to use for this worker.
    Returns:
        None; saves results to disk.
    """
    print(f"Processing region: {region}")

    # Initialize DESeq2 inference method
    inference = DefaultInference(cpus_per_worker)

    # Load adata
    adata_fp = f'{adata_input_prefix}_{region}.h5ad'
    print(f"Loading data from: {adata_fp}")
    adata = sc.read_h5ad(adata_fp)

    # Initialize DESeq2 dataset
    print("Initializing DESeq2 dataset")
    dds = initialize_dds(adata=adata, inference=inference)

    print("Running DESeq2 analysis")
    results = run_deseq2(dds=dds, inference=inference)

    # Save results
    results_fp = f'{results_prefix}_{region}.tsv'
    print(f"Saving results to: {results_fp}")
    results.to_csv(results_fp, sep='\t')

    # Identify significant genes
    significant_genes = identify_significant_genes(results)

    # Clean memory
    del adata, dds, results
    gc.collect()

    return significant_genes

def main(args: argparse.Namespace):

    # Extract arguments
    regions = args.regions.split()
    adata_input_prefix = args.adata_input_prefix
    results_prefix = args.results_prefix

    # Number of parallel workers
    n_workers = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))

    # Allocate CPUs
    total_cpus = os.cpu_count()
    cpus_per_worker = max(1, total_cpus // n_workers)
    
    print(f"Total CPUs: {total_cpus}, CPUs per worker: {cpus_per_worker}")

    # Store significant genes DataFrames with region info
    significant_genes_by_region = {}

    # Process each region in parallel
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {}
        for region in regions:
            future = executor.submit(
                process_region,
                region=region, adata_input_prefix=adata_input_prefix, results_prefix=results_prefix, cpus_per_worker=cpus_per_worker)
            futures[future] = region

        # Wait for tasks as they complete (not in order)
        for future in as_completed(futures):
            region = futures[future]
            try:
                significant_genes = future.result()
                significant_genes_by_region[region] = significant_genes
                print(f"✓ Region {region} completed")
            except Exception as e:
                print(f"✗ Region {region} failed: {e}")
                traceback.print_exc()
            finally:
                del future

    # Concatenate all significant genes DataFrames with region column
    if significant_genes_by_region:
        dfs_to_concat = []
        for region, df in significant_genes_by_region.items():
            df = df.copy()
            df['Region'] = region
            dfs_to_concat.append(df)
        
        all_significant_genes = pd.concat(dfs_to_concat, ignore_index=False)
        
        # Save concatenated results
        output_fp = f'{results_prefix}_all_significant_genes.tsv'
        all_significant_genes.to_csv(output_fp, sep='\t')
        print(f"\nAll significant genes saved to: {output_fp}")
        print(f"Total significant genes across all regions: {len(all_significant_genes)}")
    else:
        print("No significant genes found in any region")
    
    # Clean up memory
    del significant_genes_by_region
    gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DESeq2 analysis on multiple regions in parallel.")
    parser.add_argument("--regions", type=str, required=True, help="Path to the regions file.")
    parser.add_argument(
        "--adata-input-prefix",
        type=str, required=True,
        help="Prefix for input AnnData files."
    )
    parser.add_argument(
        "--results-prefix",
        type=str,
        required=True,
        help="Prefix for output DESeq2 results files."
    )

    args = parser.parse_args()
    main(args)
